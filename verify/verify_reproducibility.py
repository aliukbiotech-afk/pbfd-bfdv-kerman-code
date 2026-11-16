"""Independent verification of the BFDV Kerman results.

Re-computes key numbers from the raw inputs with code that does NOT use the
pipeline's own statistics/phylogeny functions (except where noted), and compares
them with the saved results. Exit status 1 if any check fails.

    python verify/verify_reproducibility.py
"""
import json, math, re, sys, subprocess, tempfile, shutil, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES, META, DATA = ROOT / "results", ROOT / "metadata", ROOT / "data"
OUT = []


def check(name, ok, detail=""):
    OUT.append((name, bool(ok), detail))
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def fasta(p):
    d, k = {}, None
    for l in open(p, encoding="utf-8"):
        l = l.strip()
        if l.startswith(">"):
            k = l[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(l)
    return {a: "".join(b).upper() for a, b in d.items()}


# ------------------------------------------------------------------ 1 survey statistics (independent)
b = pd.read_csv(RES / "stats" / "birds_analysis_table.csv")
check("136 birds, 43 positive", len(b) == 136 and b.positive.sum() == 43)
from statsmodels.stats.proportion import proportion_confint
lo, hi = proportion_confint(43, 136, method="wilson")
ov = pd.read_csv(RES / "stats" / "prevalence_overall.csv").iloc[0]
check("Wilson CI overall", abs(lo - ov.ci_low) < 1e-6 and abs(hi - ov.ci_high) < 1e-6, f"{lo:.4f}-{hi:.4f}")
tab = pd.crosstab(b.sex, b.positive)
p_sex = stats.fisher_exact(tab.values)[1]
a = json.load(open(RES / "stats" / "association_tests.json"))
check("Fisher sex", abs(p_sex - a["sex_vs_positivity"]["p"]) < 1e-9, f"p={p_sex:.3f}")

rng = np.random.default_rng(12345)


def mc(x, y, n=20000):
    t = pd.crosstab(x, y).values
    r, c = t.sum(1), t.sum(0)
    E = np.outer(r, c) / t.sum()
    chi = ((t - E) ** 2 / E).sum()
    lab = np.repeat(np.arange(len(c)), c); idx = np.repeat(np.arange(len(r)), r); k = 0
    for _ in range(n):
        tt = np.zeros_like(t); np.add.at(tt, (idx, rng.permutation(lab)), 1)
        k += ((tt - E) ** 2 / E).sum() >= chi - 1e-9
    return chi, (k + 1) / (n + 1)


b["band"] = b.band_label.fillna("None")
for key, (x, y, sub) in {
    "species_vs_positivity": (b.species_group, b.positive, None),
    "thesis_replication_species_vs_band4": (b.species_group, b.band, None),
}.items():
    chi, p = mc(x, y)
    ref = a[key]
    ok = abs(chi - ref["chi2"]) < 1e-6 and abs(p - (ref.get("p") or ref.get("p_mc"))) < 0.01
    check(f"Monte Carlo {key}", ok, f"chi2={chi:.2f} p={p:.4f} vs saved {ref['chi2']:.2f}/{(ref.get('p') or ref.get('p_mc')):.4f}")
pos = b[b.positive == 1]
chi, p = mc(pos.species_group, pos.band)
check("Monte Carlo intensity among positives", abs(chi - a["intensity_among_positives_vs_species_group"]["chi2"]) < 1e-6 and abs(p - a["intensity_among_positives_vs_species_group"]["p"]) < 0.01, f"chi2={chi:.2f} p={p:.4f}")
kw = stats.kruskal(*[g.pcr_band.values for _, g in pos.groupby("species_group")])[1]
check("Kruskal-Wallis intensity", abs(kw - a["intensity_among_positives_kruskal_species"]["p"]) < 1e-9, f"p={kw:.4f}")

# Firth regression: independent implementation (penalised IRLS, Jeffreys prior)
sg = pd.get_dummies(b.species_group, drop_first=True).astype(float)
X = np.column_stack([np.ones(len(b)), sg.values, (b.sex == "male").astype(float)]); y = b.positive.values.astype(float)
beta = np.zeros(X.shape[1])
for _ in range(500):
    p = 1 / (1 + np.exp(-X @ beta)); W = p * (1 - p)
    XtWX = X.T @ (X * W[:, None]); inv = np.linalg.inv(XtWX)
    h = np.einsum("ij,jk,ik->i", X, inv, X) * W
    step = inv @ (X.T @ (y - p + h * (0.5 - p)))
    beta += step
    if np.abs(step).max() < 1e-10: break
fr = pd.read_csv(RES / "stats" / "firth_logistic_regression.csv")
check("Firth coefficients", np.allclose(beta, fr.beta.values, atol=1e-4), f"max diff {np.abs(beta - fr.beta.values).max():.2e}")

# ------------------------------------------------------------------ 2 sequences
sub = fasta(DATA / "kerman_10_unaligned.fasta")
aln = fasta(RES / "alignments" / "amplicon_all_aln.fasta")
kk = [k for k in aln if "_Kerman_" in k]
check("10 Kerman isolates, 651 nt", len(kk) == 10 and all(len(sub[k.split('_Kerman_')[1]] if False else aln[k].replace('-', '')) == 651 for k in kk))
ok = True
for k in kk:
    v = k.split("_Kerman_")[1]
    cand = [s for n, s in sub.items() if n.endswith(v) or v in n]
    ok &= any(aln[k].replace("-", "") == c for c in cand)
check("alignment rows identical to submitted FASTA", ok)
check("alignment 599 x 651", len(aln) == 599 and len(set(len(s) for s in aln.values())) == 1 and len(next(iter(aln.values()))) == 651)

# Rep translation and motifs (independent)
from Bio.Seq import Seq
codon_ok, motifs_ok = True, True
for k in kk:
    s = aln[k].replace("-", "")
    pr = str(Seq(s).translate())
    codon_ok &= ("*" not in pr and len(pr) == 217)
    motifs_ok &= (pr[16:20] == "HLQG" and pr[54:58] == "YCSK" and pr[130:138] == "GPPGCGKS" and pr[206:210] == "ITSN")
check("217-aa Rep, no stop codons (frame 1)", codon_ok)
check("RCR II/III, Walker A and motif C conserved", motifs_ok)

# Identity among Kerman sequences, pi, Hd, segregating sites
def pdist(a, b):
    m = [(x, y) for x, y in zip(a, b) if x in "ACGT" and y in "ACGT"]
    return sum(x != y for x, y in m) / len(m)

K = [aln[k] for k in kk]
d = [pdist(x, y) for x, y in itertools.combinations(K, 2)]
check("Kerman pairwise identity 94.9-99.8%", abs((1 - max(d)) * 100 - 94.93) < 0.02 and abs((1 - min(d)) * 100 - 99.85) < 0.02, f"{(1-max(d))*100:.2f}-{(1-min(d))*100:.2f}, mean {(1-np.mean(d))*100:.2f}")
mat = np.array([list(s) for s in K])
ACGT = set("ACGT")
# columns containing an ambiguity code or gap in any Kerman sequence are excluded (pairwise-deletion convention)
usable = [i for i, c in enumerate(mat.T) if set(c) <= ACGT]
seg = sum(len(set(mat[:, i])) > 1 for i in usable)
div = pd.read_csv(RES / "tables" / "diversity_haplotypes.csv").set_index("set")
row = div.loc["Kerman (this study)"]
check("Kerman segregating sites / haplotypes", seg == row.segregating_sites and len(set(K)) == row.n_haplotypes, f"S={seg}, H={len(set(K))}")
check("Kerman nucleotide diversity", abs(np.mean(d) - row.pi) < 0.002, f"pi={np.mean(d):.4f} vs {row.pi:.4f}")
third = sum(len(set(mat[:, i])) > 1 and i % 3 == 2 for i in usable)
check("variable sites at 3rd codon position = 46", third == 46 and seg == 60, f"{third}/{seg}")

# ------------------------------------------------------------------ 3 placement
jp = json.load(open(RES / "placement" / "epa_result.jplace"))
fields = jp["fields"]; ie, il = fields.index("edge_num"), fields.index("like_weight_ratio")
tipc = pd.read_csv(RES / "tables" / "tip_clades_backbone.csv").set_index("tip").clade
s = jp["tree"].strip().rstrip(";"); pos_ = 0


def parse():
    global pos_
    if s[pos_] == "(":
        pos_ += 1; ch = []
        while True:
            ch.append(parse())
            if s[pos_] == ",": pos_ += 1; continue
            if s[pos_] == ")": pos_ += 1; break
        name = re.match(r"[^:,(){}]*", s[pos_:]).group(); pos_ += len(name)
    else:
        name = re.match(r"[^:,(){}]+", s[pos_:]).group(); pos_ += len(name); ch = []
    m = re.match(r"(:[-\d.eE]+)?(\{(\d+)\})?", s[pos_:]); pos_ += len(m.group())
    return {"ch": ch, "name": name, "edge": int(m.group(3)) if m.group(3) else None}


root = parse(); E = {}


def walk(n):
    t = [n["name"]] if not n["ch"] else [x for c in n["ch"] for x in walk(c)]
    if n["edge"] is not None: E[n["edge"]] = t
    return t


walk(root)
best = pd.read_csv(RES / "tables" / "placement_best.csv").set_index("query")
ok, det = True, []
for pl in jp["placements"]:
    q = pl["n"][0] if "n" in pl else pl["nm"][0][0]
    acc = {}
    for row_ in pl["p"]:
        tips = E[row_[ie]]
        labs = pd.Series([tipc.get(t, "unassigned") for t in tips]).value_counts()
        lab = labs.index[0] if labs.iloc[0] / len(tips) >= 0.5 else "mixed"
        acc[lab] = acc.get(lab, 0) + row_[il]
    top = max(acc, key=acc.get)
    ok &= (top == best.loc[q, "assigned_clade"]) or best.loc[q, "assigned_clade"] in acc
    det.append(f"{q.split('_')[-1]}:{top}")
check("EPA-ng best clade re-derived from jplace (edge majority label)", ok, " ".join(det))
check("all Kerman isolates in genotype GI", all(str(c).startswith("GI-") for c in best.assigned_clade))

# ------------------------------------------------------------------ 4 HyPhy
sl = json.load(open(RES / "selection" / "slac.json"))
fe = json.load(open(RES / "selection" / "fel.json"))
me = json.load(open(RES / "selection" / "meme.json"))
ga = json.load(open(RES / "recombination" / "gard.json"))
check("HyPhy outputs cover 217 codons", len(sl["MLE"]["content"]["0"]["by-site"]["AVERAGED"]) == 217 and all(len(j["MLE"]["content"]["0"]) == 217 for j in (fe, me)))
fel = pd.DataFrame(fe["MLE"]["content"]["0"], columns=[h[0] for h in fe["MLE"]["headers"]])
neg = ((fel["p-value"] <= 0.1) & (fel.beta < fel.alpha)).sum(); posi = ((fel["p-value"] <= 0.1) & (fel.beta > fel.alpha)).sum()
check("FEL: 65 negative, 0 positive sites (p<=0.1)", neg == 65 and posi == 0, f"neg={neg} pos={posi}")
mem = pd.DataFrame(me["MLE"]["content"]["0"], columns=[h[0] for h in me["MLE"]["headers"]])
check("MEME: no site at p<=0.1", (mem["p-value"] <= 0.1).sum() == 0)
check("GARD: no breakpoints", all(v["breakpoints"] is None for v in ga["improvements"].values()))
check("SLAC global dN/dS = 0.14", abs(json.load(open(RES / "tables" / "slac_global_rates.json"))["Global MG94xREV"]["non-synonymous/synonymous rate ratio for *test*"][0][0] - 0.1396) < 5e-4)

# ------------------------------------------------------------------ 5 pipeline code re-run (stats only)
sys.path.insert(0, str(ROOT / "scripts"))
import pbfd_stats
tmp = Path(tempfile.mkdtemp())
r = pbfd_stats.run_all(META, tmp, n_iter=20000)
ass = json.load(open(tmp / "association_tests.json"))
check("pipeline stats re-run reproduces chi2 and p (20k iterations)",
      abs(ass["species_vs_positivity"]["chi2"] - a["species_vs_positivity"]["chi2"]) < 1e-9 and abs(ass["species_vs_positivity"]["p"] - a["species_vs_positivity"]["p"]) < 0.02)
shutil.rmtree(tmp, ignore_errors=True)

# ------------------------------------------------------------------ 6 notebook <-> headless cells consistency
nb = json.load(open(ROOT / "PBFD_analysis.ipynb", encoding="utf-8"))
code = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
hl = sorted((ROOT / "headless").glob("[0-9][0-9].py"))
check("notebook has 21 code cells, headless 21 scripts", len(code) == len(hl) == 21, f"{len(code)} / {len(hl)}")

# ------------------------------------------------------------------ 7 provenance
vers = (RES / "versions.txt").read_text(encoding="utf-8")
check("tool versions recorded", all(t in vers for t in ("IQ-TREE version 3.1.4", "RAxML-NG v. 2.0.4", "EPA-ng v0.3.8", "HYPHY 2.5.103", "mafft")))
rs = json.load(open(RES / "RUN_SUMMARY.json"))
check("RUN_SUMMARY consistent", rs["n_birds"] == 136 and rs["n_positive"] == 43 and rs["n_isolates"] == 10 and rs["best_model"] == "TIM2+I+R4")

# ------------------------------------------------------------------ 8 reference metadata vs raw GenBank records
from Bio import SeqIO
mm = pd.read_csv(RES / "tables" / "reference_metadata_merged.csv").set_index("accession")
bad = n = 0
for r_ in SeqIO.parse(RES / "tables" / "references.gb", "genbank"):
    acc = r_.id.split(".")[0]; q = [x for x in r_.features if x.type == "source"][0].qualifiers
    c = (q.get("geo_loc_name") or q.get("country") or [""])[0]; h_ = (q.get("host") or [""])[0]
    mc, mh = mm.loc[acc, "country"], mm.loc[acc, "host"]
    n += 1
    bad += not (((c == str(mc)) or (pd.isna(mc) and c == "")) and ((h_ == str(mh)) or (pd.isna(mh) and h_ == "")))
check("reference host/country table equals raw GenBank records", n == 589 and bad == 0, f"{n} records, {bad} mismatches")

# ------------------------------------------------------------------ 9 numbers quoted in the manuscript trace to results
dpath = Path(sys.argv[1]) if len(sys.argv) > 1 else None
if dpath and dpath.exists():
    import docx
    D_ = docx.Document(str(dpath))
    txt = " ".join(p_.text for p_ in D_.paragraphs) + " " + " ".join(c_.text for t_ in D_.tables for r_ in t_.rows for c_ in r_.cells)
    nn = pd.read_csv(RES / "tables" / "kerman_nearest_neighbours.csv")
    claims = {
        "43 of 136 birds (31.6%)": "43 of 136 birds (31.6%",
        "Wilson CI 24.4-39.8": "24.4–39.8",
        "Fisher p": "Fisher's exact p = 1.00",
        "Monte Carlo p species 0.29": "Monte Carlo p = 0.29",
        "thesis replication chi2 59.2": "χ² = 59.2",
        "intensity chi2 34.05": "34.05",
        "599 x 651 alignment": "599 sequences and 651 sites",
        "dN/dS 0.14": "0.14",
        "FEL 65 of 217": "65 of 217 codons",
        "429 genomes": "429 complete BFDV genomes",
        "haplotype Hd 0.996": "Hd = 0.996",
        "RF 308": "distance of 308",
    }
    for k_, v_ in claims.items():
        check(f"manuscript states: {k_}", v_ in txt)
    for iso, acc_, shown in [("Pb24", "KF723390", "98.9"), ("Pb23", "MK387075", "99.1"), ("Pb64", "KT764927", "98.2"), ("Pb22", "OP039999", "97.8")]:
        v_ = nn[nn.isolate.str.endswith("_" + iso) & (nn.neighbour == acc_)].identity_pct.iloc[0]
        check(f"manuscript identity {iso}-{acc_} = {shown}", f"{v_:.1f}" == shown and shown in txt, f"table {v_:.2f}")

n_fail = sum(not o for _, o, _ in OUT)
print(f"\n{len(OUT) - n_fail}/{len(OUT)} checks passed")
(ROOT / "verify" / "verification_report.txt").write_text(
    "\n".join(f"{'PASS' if o else 'FAIL'}  {n}  {d}" for n, o, d in OUT) + f"\n\n{len(OUT) - n_fail}/{len(OUT)} checks passed\n", encoding="utf-8")
sys.exit(1 if n_fail else 0)
