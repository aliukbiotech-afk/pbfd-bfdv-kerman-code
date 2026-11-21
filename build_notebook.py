"""Generate PBFD_analysis.ipynb (Colab).

    python build_notebook.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
cells = []
md = lambda s: cells.append(dict(cell_type="markdown", metadata={}, source=s.strip("\n")))
code = lambda s: cells.append(dict(cell_type="code", metadata={}, execution_count=None, outputs=[], source=s.strip("\n")))

md(r"""
# PBFD – Kerman: reproducible re-analysis (Colab)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aliukbiotech-afk/pbfd-bfdv-kerman-code/blob/main/PBFD_analysis.ipynb)

Subclinical beak and feather disease virus (BFDV) in 136 pet psittacines, Kerman, Iran: 10 partial ORF1/Rep sequences (PX634653–PX634662).

**How to run:** click the *Open in Colab* badge, enter an e-mail address for NCBI Entrez in the config cell (NCBI requires one), then `Runtime → Run all`. The inputs (`data/`, `metadata/`, `scripts/`) are cloned automatically from https://github.com/aliukbiotech-afk/pbfd-bfdv-kerman-code; nothing has to be uploaded. At the end the notebook downloads **`results.zip`**.

Runtime is about 1.5–3 h on a standard CPU runtime. Set `FAST = True` in the config cell for a 20–30 min test run with reduced bootstraps.

| Section | Content |
|---|---|
| 0 | Setup: bioconda tools, inputs, config |
| 1 | Sanger chromatogram QC |
| 2 | Reference genomes and metadata from GenBank |
| 3 | Amplicon extraction and alignment |
| 4 | ML tree (ModelFinder, UFBoot + SH-aLRT) and thesis-model replication |
| 5 | Full-genome backbone, Shah 2023 sub-clades, EPA-ng placement of the partial sequences |
| 6 | Distances, nearest neighbours, haplotypes and network |
| 7 | Rep motifs, recombination (GARD), selection (SLAC/FEL/MEME) |
| 8 | Phylogeny–trait association (host, country) |
| 9 | Epidemiological statistics |
| 10 | Figures and export |

Host and species names are read from the `FINAL_*` columns of `metadata/*.csv`; edit them there and `scripts/make_figures.py` re-renders every figure locally. Figure colours follow `FIGURE_COLOUR_KEY.md`.
""")

md("## 0. Setup")
code(r"""
#@title Config
FAST = False          #@param {type:"boolean"}
EMAIL = ""            #@param {type:"string"}
SEED = 20261103
import os, sys, pathlib, time, json, shutil, subprocess
IN_COLAB = "google.colab" in sys.modules
WORK = pathlib.Path("/content/pbfd") if IN_COLAB else pathlib.Path.cwd() / "colab_work"
WORK.mkdir(parents=True, exist_ok=True); os.chdir(WORK)
THREADS = os.cpu_count() or 2
UFBOOT = 1000; ALRT = 1000; NPERM = 1000; NMC = 100_000
if FAST:
    UFBOOT, ALRT, NPERM, NMC = 1000, 0, 200, 20_000
print("work dir:", WORK, "| threads:", THREADS, "| FAST:", FAST)
""")
code(r"""
#@title Install tools (micromamba + bioconda; ~5 min)
ENV = pathlib.Path("/content/env") if IN_COLAB else pathlib.Path.home() / "pbfd_env"
BIN = ENV / "bin"
if not (BIN / "iqtree2").exists() and not (BIN / "iqtree").exists():
    subprocess.run("curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba", shell=True, check=True)
    subprocess.run(f"micromamba create -y -q -p {ENV} -c conda-forge -c bioconda "
                   "mafft iqtree raxml-ng epa-ng hyphy", shell=True, check=True)
subprocess.run("pip -q install biopython networkx", shell=True, check=True)
IQTREE = str(BIN / "iqtree2") if (BIN / "iqtree2").exists() else str(BIN / "iqtree")
MAFFT, RAXML, EPA, HYPHY = (str(BIN / x) for x in ("mafft", "raxml-ng", "epa-ng", "hyphy"))
""")
code(r"""
#@title Get inputs (cloned from GitHub)
REPO_URL = "https://github.com/aliukbiotech-afk/pbfd-bfdv-kerman-code"
INP = WORK / "inputs"
if not IN_COLAB and (WORK.parent / "scripts" / "pbfd_phylo.py").exists():
    INP = WORK.parent                      # notebook run inside a local clone of the repository
elif not (INP / "scripts").exists():
    subprocess.run(["git", "clone", "-q", "--depth", "1", REPO_URL, str(INP)], check=True)
DATA, META = INP / "data", INP / "metadata"
sys.path.insert(0, str(INP / "scripts"))
import pbfd_phylo as P, pbfd_stats as S, pbfd_figures as F, make_figures as MF
from pbfd_names import pending_names
import numpy as np, pandas as pd
RES = WORK / "results"
for d in ["tables", "trees", "alignments", "figures", "placement", "selection", "recombination", "logs", "stats"]:
    (RES / d).mkdir(parents=True, exist_ok=True)
LOG = RES / "logs" / "commands.log"
run = lambda cmd, **kw: P.run(cmd, LOG, **kw)
print(sorted(p.name for p in DATA.iterdir())); print(pending_names(META))
""")
code(r"""
#@title Software versions
vers = {}
for name, cmd in [("mafft", f"{MAFFT} --version"), ("iqtree", f"{IQTREE} --version"), ("raxml-ng", f"{RAXML} --version"),
                  ("epa-ng", f"{EPA} --version"), ("hyphy", f"{HYPHY} --version")]:
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    vers[name] = (out.stdout + out.stderr).strip().splitlines()[0:3]
import Bio, scipy, statsmodels, matplotlib, networkx
vers.update(python=sys.version.split()[0], biopython=Bio.__version__, numpy=np.__version__, pandas=pd.__version__,
            scipy=scipy.__version__, statsmodels=statsmodels.__version__, matplotlib=matplotlib.__version__, networkx=networkx.__version__,
            seed=SEED, FAST=FAST, ufboot=UFBOOT, alrt=ALRT, permutations=NPERM, monte_carlo=NMC)
json.dump(vers, open(RES / "versions.json", "w"), indent=2)
(RES / "versions.txt").write_text("\n".join(f"{k}: {v}" for k, v in vers.items()))
print(open(RES / "versions.txt").read())
""")

md("## 1. Sanger chromatogram QC")
code(r"""
from Bio import SeqIO, Align
kerman = P.read_fasta(DATA / "kerman_10_unaligned.fasta")
hosts = pd.read_csv(META / "host_mapping.csv", dtype=str, keep_default_na=False)
al = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)
rows = []
for ab in sorted((DATA / "ab1").glob("*.ab1")):
    rec = SeqIO.read(ab, "abi"); q = np.array(rec.letter_annotations["phred_quality"])
    tr = SeqIO.read(ab, "abi-trim")
    stem = ab.stem; no = int(stem.split("-")[0])
    h = hosts[hosts.seq_no == str(no)].iloc[0]
    sub = next(v for k, v in kerman.items() if k.startswith(h.accession))
    best = max((al.align(str(x), sub)[0] for x in (tr.seq, tr.seq.reverse_complement())), key=lambda a: a.score)
    cnt = best.counts()
    cov = sum(b1 - b0 for b0, b1 in best.aligned[1]) / len(sub)
    rows.append(dict(seq_no=no, voucher=h.voucher, accession=h.accession, raw_length=len(rec), trimmed_length=len(tr),
                     mean_phred=q.mean(), q20_fraction=(q >= 20).mean(), q30_fraction=(q >= 30).mean(),
                     n_ambiguous_calls=sum(c not in "ACGT" for c in str(rec.seq)),
                     submitted_length=len(sub), submitted_covered_by_read=cov,
                     identity_read_vs_submitted=cnt.identities / max(1, cnt.identities + cnt.mismatches),
                     gaps_read_vs_submitted=cnt.gaps))
qc = pd.DataFrame(rows).sort_values("seq_no"); qc.to_csv(RES / "tables" / "sanger_qc.csv", index=False)
raw = P.read_fasta(DATA / "kerman_10_raw_basecalls.fasta")
qc["lowercase_edited_bases_in_raw"] = [sum(c.islower() for c in next((v for k, v in raw.items() if k.startswith(f"{n}-")), "")) for n in qc.seq_no]
qc.to_csv(RES / "tables" / "sanger_qc.csv", index=False); qc.round(3)
""")

md("## 2. Reference genomes and metadata (GenBank)")
code(r"""
thesis_aln = P.read_fasta(DATA / "thesis_alignment_598.fasta")
REFS = [k for k in thesis_aln if "_Kerman_" not in k]
print(len(REFS), "reference accessions used in the thesis tree")
assert EMAIL and "@" in EMAIL, "Enter your e-mail in the Config cell (NCBI requirement)"
gbmd = P.entrez_fetch(REFS, EMAIL, RES / "tables")
old = pd.read_csv(META / "reference_metadata.csv")
mdm = gbmd.merge(old[["accession", "isolation_source"]], on="accession", how="left")
shah = pd.read_csv(META / "shah2023_clades.csv")
mdm = mdm.merge(shah[["accession", "clade"]].rename(columns={"clade": "shah2023_clade"}), on="accession", how="left")
mdm.to_csv(RES / "tables" / "reference_metadata_merged.csv", index=False)
print("records:", len(mdm), "| complete genomes:", mdm.complete_genome.sum(), "| with collection year:", mdm.year.notna().sum(),
      "| with Shah-2023 clade:", mdm.shah2023_clade.notna().sum())
print(mdm.country_simple.value_counts().head(15))
""")

md("## 3. Amplicon extraction and alignment")
code(r"""
genomes = P.read_fasta(RES / "tables" / "reference_genomes_raw.fasta")
probe = max(kerman.values(), key=len)
oriented, orient_log = P.orient_and_rotate(genomes, probe)
orient_log.to_csv(RES / "tables" / "genome_orientation_log.csv", index=False)
P.write_fasta(oriented, RES / "alignments" / "reference_genomes_oriented.fasta")
amp, amp_log = P.extract_amplicon(oriented, probe)
amp_log.to_csv(RES / "tables" / "amplicon_extraction_log.csv", index=False)
bad = amp_log[(amp_log.probe_cov < 0.9)].accession.tolist()
print("genomes:", len(oriented), "| nonamer found:", orient_log.nonamer_found.sum(), "| amplicons with <90% probe coverage (dropped):", bad)
amp = {k: v for k, v in amp.items() if k not in bad}
allamp = dict(amp); allamp.update(kerman)
P.write_fasta(allamp, RES / "alignments" / "amplicon_all_unaligned.fasta")
run(f"{MAFFT} --auto --thread {THREADS} --quiet {RES/'alignments'/'amplicon_all_unaligned.fasta'} > {RES/'alignments'/'amplicon_raw_aln.fasta'}")
aln = P.trim_to_query_span(P.read_fasta(RES / "alignments" / "amplicon_raw_aln.fasta"), list(kerman))
aln = P.drop_allgap_columns(aln)
P.write_fasta(aln, RES / "alignments" / "amplicon_all_aln.fasta")
KER = list(kerman)
print("alignment:", len(aln), "sequences x", len(next(iter(aln.values()))), "columns (amplicon window spanned by all Kerman isolates)")
""")

md("## 4. Maximum-likelihood tree of the amplicon")
code(r"""
A = RES / "alignments" / "amplicon_all_aln.fasta"; T = RES / "trees"
alrt = f"-alrt {ALRT}" if ALRT else ""
run(f"{IQTREE} -s {A} -m MFP -B {UFBOOT} {alrt} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'amplicon_MFP'} -redo")
run(f"{IQTREE} -s {A} -m TIM3+F+I+G4 -B {UFBOOT} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'amplicon_TIM3'} -redo")
best_model = [l for l in open(T / "amplicon_MFP.iqtree") if l.startswith("Best-fit model")][0].split(":")[1].strip()
run(f"{IQTREE} -rf {T/'amplicon_MFP.treefile'} {T/'amplicon_TIM3.treefile'} --prefix {T/'rf_MFP_vs_TIM3'} -redo")
print("best-fit model (BIC):", best_model)
print(open(T / "rf_MFP_vs_TIM3.rfdist").read())
""")

md("## 5. Full-genome backbone, Shah 2023 sub-clades and placement")
code(r"""
full = {k: v for k, v in oriented.items() if 1900 <= len(v) <= 2100}
P.write_fasta(full, RES / "alignments" / "backbone_genomes_unaligned.fasta")
run(f"{MAFFT} --auto --thread {THREADS} --quiet {RES/'alignments'/'backbone_genomes_unaligned.fasta'} > {RES/'alignments'/'backbone_genomes_aln.fasta'}")
run(f"{IQTREE} -s {RES/'alignments'/'backbone_genomes_aln.fasta'} -m MFP -B {UFBOOT} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'backbone_fullgenome_MFP'} -redo")
shah_lab = dict(zip(shah.accession, shah.clade))
bb = P.midpoint(P.read_tree(T / "backbone_fullgenome_MFP.treefile"))
bclades = P.propagate_clades(bb, shah_lab)
bclades.to_csv(RES / "tables" / "tip_clades_backbone.csv", index=False)
json.dump(bclades.attrs["clade_nodes"], open(RES / "tables" / "backbone_clade_nodes.json", "w"), indent=2)
print("backbone genomes:", len(full)); print(bclades.groupby(["clade", "source"]).size().unstack(fill_value=0))
""")
code(r"""
#@title Clade labels on the amplicon tree (same propagation, for comparison)
at = P.midpoint(P.read_tree(T / "amplicon_MFP.treefile"))
lab_all = dict(zip(bclades.tip, bclades.clade))
lab_ref = {k: v for k, v in lab_all.items() if v != "unassigned"}
aclades = P.propagate_clades(at, lab_ref)
aclades.to_csv(RES / "tables" / "tip_clades_amplicon.csv", index=False)
aclades[aclades.tip.isin(KER)]
""")
code(r"""
#@title EPA-ng placement of the 10 partial sequences on the full-genome backbone
PL = RES / "placement"
refs_bb = [t.name for t in bb.get_terminals() if t.name in aln]
P.write_fasta({k: aln[k] for k in refs_bb}, PL / "ref_amplicon_msa.fasta")
P.write_fasta({k: aln[k] for k in KER}, PL / "query.fasta")
import re
nwk = open(T / "backbone_fullgenome_MFP.treefile").read()
nwk = re.sub(r"\)[0-9.\/]+:", "):", nwk)                       # strip support labels
from Bio import Phylo
tt = Phylo.read(__import__("io").StringIO(nwk), "newick")
for tip in [x for x in tt.get_terminals() if x.name not in set(refs_bb)]:
    tt.prune(tip)
Phylo.write(tt, PL / "backbone_pruned.nwk", "newick")
(PL / "backbone_pruned.nwk").write_text(re.sub(r"\)[0-9.]+:", "):", (PL / "backbone_pruned.nwk").read_text()))
run(f"{RAXML} --evaluate --msa {PL/'ref_amplicon_msa.fasta'} --tree {PL/'backbone_pruned.nwk'} --model GTR+G --prefix {PL/'eval'} --threads {THREADS} --redo --seed {SEED}")
run(f"{EPA} --ref-msa {PL/'ref_amplicon_msa.fasta'} --tree {PL/'eval.raxml.bestTree'} --query {PL/'query.fasta'} "
    f"--model {PL/'eval.raxml.bestModel'} --outdir {PL} --redo --threads {THREADS}")
ptree, pl = P.parse_jplace(PL / "epa_result.jplace")
tipc = dict(zip(bclades.tip, bclades.clade))
pl, summ, best = P.placement_clade_summary(ptree, pl, tipc)
pl.to_csv(RES / "tables" / "placement_all.csv", index=False)
summ.to_csv(RES / "tables" / "placement_lwr_by_clade.csv", index=False)
best.to_csv(RES / "tables" / "placement_best.csv", index=False)
best
""")
code(r"""
#@title Cross-check: IQ-TREE with the backbone as topological constraint
run(f"{IQTREE} -s {A} -g {PL/'backbone_pruned.nwk'} -m {best_model} -B {UFBOOT} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'amplicon_constrained'} -redo")
ct = P.midpoint(P.read_tree(T / "amplicon_constrained.treefile"))
cc = P.propagate_clades(ct, lab_ref)
cmp = (aclades[aclades.tip.isin(KER)][["tip", "clade"]].rename(columns={"clade": "amplicon_tree"})
       .merge(cc[cc.tip.isin(KER)][["tip", "clade"]].rename(columns={"clade": "constrained_tree"}), on="tip")
       .merge(best[["query", "assigned_clade", "clade_LWR"]].rename(columns={"query": "tip", "assigned_clade": "EPA_ng"}), on="tip"))
cmp["agree_all"] = (cmp.amplicon_tree == cmp.constrained_tree) & (cmp.constrained_tree == cmp.EPA_ng)
cmp.to_csv(RES / "tables" / "kerman_clade_assignment_comparison.csv", index=False); cmp
""")

md("## 6. Distances, nearest neighbours, diversity, haplotype network")
code(r"""
D = P.p_distance_matrix(aln, KER)
mdd = mdm.set_index("accession")
rows = []
for k in KER:
    for acc, d in D.loc[k].drop(KER).nsmallest(10).items():
        r = mdd.loc[acc] if acc in mdd.index else {}
        rows.append(dict(isolate=k, neighbour=acc, p_distance=d, identity_pct=100 * (1 - d), country=r.get("country_simple", ""),
                         host=r.get("host", ""), year=r.get("year_or_submit", ""), clade=tipc.get(acc, "")))
near = pd.DataFrame(rows); near.to_csv(RES / "tables" / "kerman_nearest_neighbours.csv", index=False)
# mean distance of each isolate to each sub-clade
dc = D.drop(columns=KER).T.assign(clade=lambda x: x.index.map(lambda a: tipc.get(a, "unassigned"))).groupby("clade").mean().T
(100 * (1 - dc)).round(2).to_csv(RES / "tables" / "kerman_mean_identity_by_clade.csv")
IRAN = [a for a in aln if a in mdd.index and mdd.loc[a, "country_simple"] == "Iran"]
Dki = P.p_distance_matrix({k: aln[k] for k in KER + IRAN})
(100 * (1 - Dki)).round(2).to_csv(RES / "tables" / "identity_matrix_kerman_iran.csv")
print("Iranian reference strains:", len(IRAN))
near.groupby("isolate").head(3)
""")
code(r"""
div = []
for name, members in [("Kerman (this study)", KER), ("Iran (all, incl. this study)", KER + IRAN), ("All sequences", list(aln))]:
    _, st = P.haplotypes(aln, members); st.insert(0, "set", name); div.append(st)
div = pd.concat(div); div.to_csv(RES / "tables" / "diversity_haplotypes.csv", index=False)
# network: Kerman isolates + their nearest reference each
net = list(dict.fromkeys(KER + near.groupby("isolate").head(1).neighbour.tolist()))
mem = pd.DataFrame(dict(tip=net, group=["this study (Kerman)" if t in KER else ("Iran (other)" if t in IRAN else "nearest non-Iranian") for t in net]))
hap, _ = P.haplotypes(aln, net); mem["haplotype"] = mem.tip.map(hap["hap_of"])
mem.to_csv(RES / "tables" / "network_members.csv", index=False)
# focus set for Fig 3: Kerman + Iranian + 3 nearest per isolate + 2 labelled representatives per sub-clade
reps = []
for c, g in aclades[aclades.source == "Shah2023"].groupby("clade"):
    reps += g.tip.head(2).tolist()
focus = list(dict.fromkeys(KER + IRAN + near.groupby("isolate").head(3).neighbour.tolist() + reps))
(RES / "tables" / "focus_tips.txt").write_text("\n".join(focus))
div
""")

md("## 7. Rep motifs, recombination and selection")
code(r"""
mot = P.motif_table({k: aln[k].replace("-", "") for k in KER})
mot.to_csv(RES / "tables" / "rep_motifs_kerman.csv", index=False)
# subset for recombination / selection: Kerman + Iranian + 5 nearest per isolate (in-frame, no stops)
sub_names = list(dict.fromkeys(KER + IRAN + near.groupby("isolate").head(5).neighbour.tolist()))
subnt = {k: aln[k].replace("-", "") for k in sub_names}
subnt = {k: v for k, v in subnt.items() if "*" not in P.translate(v)[:-1] and len(v) >= 600}
print("selection/recombination subset:", len(subnt), "sequences")
SEL = RES / "selection"
cod = P.codon_alignment(subnt, MAFFT, SEL, LOG)
cod = {k: v[: len(v) // 3 * 3] for k, v in cod.items()}
P.write_fasta(cod, SEL / "subset_codon_aln.fasta")
P.write_fasta(cod, RES / "recombination" / "subset_codon_aln_for_RDP5.fasta")   # for optional RDP5 run on Windows
run(f"{IQTREE} -s {SEL/'subset_codon_aln.fasta'} -m GTR+F+G4 -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {SEL/'subset'} -redo")
# motif codon ranges in the codon alignment (for Fig S3)
ref_aa = P.translate(cod[KER[0]].replace("-", ""))
import re as _re
cols = [i for i in range(0, len(cod[KER[0]]), 3) if cod[KER[0]][i:i+3] != "---"]
mrows = []
for k, pat in P.REP_MOTIFS.items():
    m = _re.search(pat, ref_aa)
    if m:
        mrows.append(dict(motif=k, start=cols[m.start()] // 3 + 1, end=cols[m.end() - 1] // 3 + 1))
pd.DataFrame(mrows).to_csv(RES / "tables" / "rep_motif_codons.csv", index=False)
mot
""")
code(r"""
#@title GARD recombination screen (HyPhy)
REC = RES / "recombination"
try:
    run(f"{HYPHY} gard --alignment {SEL/'subset_codon_aln.fasta'} --type nucleotide --rv GDD --rate-classes 3 "
        f"--output {REC/'gard.json'} --output-lf {REC/'gard_best.nex'} CPU={THREADS}")
    G = json.load(open(REC / "gard.json"))
    bp = G.get("breakpointData", {})
    print("GARD: potential breakpoints:", G.get("potentialBreakpoints"), "| best model breakpoints:", len(bp) - 1 if bp else 0)
    print(json.dumps({k: G[k] for k in G if k in ("improvements", "siteBreakPointSupport")}, indent=1)[:3000])
except Exception as e:
    print("GARD failed:", e)
""")
code(r"""
#@title Selection: SLAC (global dN/dS), FEL and MEME (per site)
tree_sel = SEL / "subset.treefile"
for m in ["slac", "fel", "meme"]:
    try:
        run(f"{HYPHY} {m} --alignment {SEL/'subset_codon_aln.fasta'} --tree {tree_sel} --output {SEL/(m+'.json')} CPU={THREADS}")
    except Exception as e:
        print(m, "failed:", e)
sites = []
for m in ["fel", "meme"]:
    if (SEL / f"{m}.json").exists():
        sites.append(P.hyphy_json_sites(SEL / f"{m}.json", m.upper()))
sites = pd.concat(sites) if sites else pd.DataFrame()
sites.to_csv(RES / "tables" / "selection_sites.csv", index=False)
if (SEL / "slac.json").exists():
    J = json.load(open(SEL / "slac.json"))
    fits = J.get("fits", {})
    glob = {k: v.get("Rate Distributions") for k, v in fits.items() if "Rate Distributions" in v}
    json.dump(glob, open(RES / "tables" / "slac_global_rates.json", "w"), indent=1, default=str); print(json.dumps(glob, indent=1, default=str)[:1500])
for meth in ["FEL", "MEME"]:
    s = sites[sites.method == meth] if len(sites) else pd.DataFrame()
    pc = next((c for c in s.columns if c.lower().startswith("p-value")), None)
    if pc:
        sig = s[s[pc].astype(float) <= 0.1]
        print(meth, "sites p<=0.1:", len(sig))
        if meth == "FEL":
            a, b = [c for c in s.columns if c.startswith("alpha")][0], [c for c in s.columns if c.startswith("beta")][0]
            print("   purifying:", int((sig[b].astype(float) < sig[a].astype(float)).sum()), "| diversifying:", int((sig[b].astype(float) > sig[a].astype(float)).sum()))
""")

md("## 8. Phylogeny–trait association")
code(r"""
state_country = {a: ("Iran" if a in KER else mdd.loc[a, "country_simple"]) for a in aln if a in KER or a in mdd.index}
state_country = {k: v for k, v in state_country.items() if isinstance(v, str) and v.strip()}
def genus(h):
    h = str(h or "").strip()
    return h.split()[0].capitalize() if h and h.lower() not in ("not found", "nan", "") and h[0].isupper() else None
state_host = {a: genus(mdd.loc[a, "host"]) for a in aln if a in mdd.index}
state_host = {k: v for k, v in state_host.items() if v}
ta = []
for name, st in [("country", state_country), ("host_genus_refs_only", state_host)]:
    r = P.trait_association(P.midpoint(P.read_tree(T / "amplicon_MFP.treefile")), st, NPERM, mc_values=None)
    r.insert(0, "trait", name); ta.append(r)
ta = pd.concat(ta); ta.to_csv(RES / "tables" / "trait_association.csv", index=False); ta
""")

md("## 9. Epidemiological statistics")
code(r"""
st = S.run_all(META, RES / "stats", NMC)
print(st["prevalence"]["overall"].to_string(index=False))
print(st["prevalence"]["by_species"][["group", "label", "n", "positive", "prevalence", "ci_low", "ci_high"]].round(3).to_string(index=False))
for k, v in st["tests"].items():
    if isinstance(v, dict):
        print(f"- {k}:", {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "table"})
st["regression"].round(3)
""")

md("## 10. Figures and export")
code(r"""
MF.render_all(RES, META)
from IPython.display import Image, display
for f in sorted((RES / "figures").glob("Fig*.png")):
    print(f.name); display(Image(str(f), width=900))
""")
code(r"""
summary = dict(n_birds=136, n_positive=int(st["birds"].positive.sum()), n_isolates=len(KER),
               amplicon_alignment=[len(aln), len(next(iter(aln.values())))], best_model=best_model,
               placement=best[["query", "assigned_clade", "clade_LWR"]].to_dict("records"),
               pending_names=pending_names(META))
json.dump(summary, open(RES / "RUN_SUMMARY.json", "w"), indent=2, default=str)
shutil.make_archive(str(WORK / "results"), "zip", RES)
print(json.dumps(summary, indent=1, default=str)[:3000])
if IN_COLAB:
    from google.colab import files; files.download(str(WORK / "results.zip"))
""")

nb = dict(nbformat=4, nbformat_minor=5, cells=cells,
          metadata=dict(kernelspec=dict(name="python3", display_name="Python 3"), language_info=dict(name="python"),
                        colab=dict(provenance=[], toc_visible=True)))
for c in nb["cells"]:
    c["source"] = c["source"].splitlines(keepends=True)
json.dump(nb, open(HERE / "PBFD_analysis.ipynb", "w", encoding="utf-8"), indent=1, ensure_ascii=False)

print("wrote", HERE / "PBFD_analysis.ipynb")
