"""Sequence / phylogenetic helpers for the PBFD Kerman analysis (used by the Colab notebook).

External programs (installed in the notebook via micromamba/bioconda):
mafft, iqtree2, raxml-ng, epa-ng, hyphy, macse (optional).
"""
from __future__ import annotations
import io, json, re, subprocess, shlex, time, itertools, math, collections
from pathlib import Path
import numpy as np
import pandas as pd

SEED = 20261103
NONAMER = "TAGTATTAC"          # BFDV origin-of-replication nonanucleotide (virion strand)


# ------------------------------------------------------------------ process helpers
def run(cmd: str, log: Path | None = None, cwd: Path | None = None, check: bool = True) -> str:
    t0 = time.time()
    p = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    out = p.stdout + ("\n[stderr]\n" + p.stderr if p.stderr.strip() else "")
    if log:
        with open(log, "a", encoding="utf-8") as fh:
            fh.write(f"\n$ {cmd}\n# exit={p.returncode} time={time.time() - t0:.1f}s\n{out[-20000:]}\n")
    if check and p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd}\n{out[-3000:]}")
    return out


def read_fasta(path) -> dict[str, str]:
    seqs, name = collections.OrderedDict(), None
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return collections.OrderedDict((k, "".join(v)) for k, v in seqs.items())


def write_fasta(seqs: dict, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn-", "TGCAYRMKVHDBNtgcayrmkvhdbn-"))[::-1]


# ------------------------------------------------------------------ NCBI
def entrez_fetch(accessions, email: str, out_dir: Path, batch: int = 100) -> pd.DataFrame:
    from Bio import Entrez, SeqIO
    Entrez.email = email or None          # blank -> Biopython warns but still works
    Entrez.tool = "pbfd_kerman_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    gb_path = out_dir / "references.gb"
    if not gb_path.exists():
        with open(gb_path, "w") as fh:
            for i in range(0, len(accessions), batch):
                ids = ",".join(accessions[i:i + batch])
                for attempt in range(4):
                    try:
                        h = Entrez.efetch(db="nuccore", id=ids, rettype="gb", retmode="text")
                        fh.write(h.read()); h.close(); break
                    except Exception as e:  # network hiccup
                        time.sleep(3 * (attempt + 1))
                else:
                    raise RuntimeError(f"Entrez failed for batch starting {accessions[i]}")
                time.sleep(0.4)
    rows, genomes = [], collections.OrderedDict()
    for rec in SeqIO.parse(gb_path, "genbank"):
        src = next((f for f in rec.features if f.type == "source"), None)
        q = src.qualifiers if src else {}
        g = lambda k: q.get(k, [""])[0]
        date = g("collection_date")
        m = re.search(r"(19|20)\d\d", date)
        acc = rec.id.split(".")[0]
        rows.append(dict(accession=acc, version=rec.id, length=len(rec.seq), definition=rec.description,
                         host=g("host"), country=g("geo_loc_name") or g("country"), collection_date=date,
                         year=int(m.group()) if m else np.nan, isolate=g("isolate") or g("strain"),
                         complete_genome=("complete genome" in rec.description.lower()),
                         submit_year=int(rec.annotations.get("date", "01-JAN-1900")[-4:])))
        genomes[acc] = str(rec.seq).upper()
    md = pd.DataFrame(rows)
    md["country_simple"] = md.country.fillna("").str.split(":").str[0].str.strip()
    md["year_or_submit"] = md.year.fillna(md.submit_year)
    md.to_csv(out_dir / "reference_genbank_metadata.csv", index=False)
    write_fasta(genomes, out_dir / "reference_genomes_raw.fasta")
    return md


def orient_and_rotate(genomes: dict, probe: str) -> tuple[dict, pd.DataFrame]:
    """Put every circular genome in the same orientation (probe = Rep amplicon) and
    rotate it to start at the origin nonanucleotide."""
    from Bio import Align
    al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                               open_gap_score=-4, extend_gap_score=-1)
    out, log = collections.OrderedDict(), []
    for acc, s in genomes.items():
        s = re.sub(r"[^ACGTRYKMSWBDHVN]", "", s.upper())
        circ = len(s) >= 1800
        f = al.score(s + s[:800] if circ else s, probe); r = al.score(revcomp(s) + revcomp(s)[:800] if circ else revcomp(s), probe)
        if r > f:
            s = revcomp(s)
        if not circ:
            out[acc] = s; log.append(dict(accession=acc, length=len(s), reversed=r > f, nonamer_found=False, score=max(f, r), circular=False))
            continue
        i = s.find(NONAMER)
        if i < 0:
            j = (s + s[:8]).find(NONAMER)
            i = j if j >= 0 else -1
        rot = s[i:] + s[:i] if i >= 0 else s
        out[acc] = rot
        log.append(dict(accession=acc, length=len(s), reversed=r > f, nonamer_found=i >= 0, score=max(f, r), circular=True))
    return out, pd.DataFrame(log)


def extract_amplicon(genomes: dict, probe: str, flank: int = 0) -> tuple[dict, pd.DataFrame]:
    """Cut the region homologous to the Kerman amplicon out of each (oriented) genome."""
    from Bio import Align
    al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                               open_gap_score=-4, extend_gap_score=-1)
    out, rows = collections.OrderedDict(), []
    for acc, s in genomes.items():
        circular = len(s) >= 1800                       # complete genomes only; partial records stay linear
        ss = s + s[:900] if circular else s
        a = al.align(ss, probe)[0]
        tb = a.aligned[0]; qb = a.aligned[1]
        t0, t1 = tb[0][0], tb[-1][1]
        q0, q1 = qb[0][0], qb[-1][1]
        # extend to cover the full probe span when local alignment stops short
        t0 = max(0, t0 - q0 - flank); t1 = min(len(ss), t1 + (len(probe) - q1) + flank)
        out[acc] = ss[t0:t1]
        rows.append(dict(accession=acc, circular=circular, start=t0, end=t1, length=t1 - t0, probe_cov=(q1 - q0) / len(probe),
                         identity_score=a.score / (2 * len(probe))))
    return out, pd.DataFrame(rows)


# ------------------------------------------------------------------ alignments / distances
def aln_array(aln: dict) -> tuple[list, np.ndarray]:
    names = list(aln)
    return names, np.array([list(aln[n].upper()) for n in names])


def p_distance_matrix(aln: dict, rows: list | None = None) -> pd.DataFrame:
    """Pairwise-deletion uncorrected p-distance (gaps/ambiguities ignored)."""
    names, A = aln_array(aln)
    valid = np.isin(A, list("ACGT"))
    idx = [names.index(r) for r in rows] if rows else range(len(names))
    D = np.zeros((len(idx), len(names)))
    for ii, i in enumerate(idx):
        both = valid[i] & valid
        diff = (A[i] != A) & both
        D[ii] = diff.sum(1) / np.maximum(both.sum(1), 1)
    return pd.DataFrame(D, index=[names[i] for i in idx], columns=names)


def trim_to_query_span(aln: dict, queries: list) -> dict:
    """Keep alignment columns spanned by every query sequence (common amplicon window)."""
    names, A = aln_array(aln)
    qi = [names.index(q) for q in queries]
    nongap = A[qi] != "-"
    cols = np.where(nongap.any(0))[0]
    lo = max(np.where(r)[0][0] for r in nongap); hi = min(np.where(r)[0][-1] for r in nongap)
    keep = [c for c in cols if lo <= c <= hi]
    return collections.OrderedDict((n, "".join(A[k, keep])) for k, n in enumerate(names))


def drop_allgap_columns(aln: dict) -> dict:
    names, A = aln_array(aln)
    keep = ~(A == "-").all(0)
    return collections.OrderedDict((n, "".join(A[k, keep])) for k, n in enumerate(names))


# ------------------------------------------------------------------ trees
def read_tree(path_or_str):
    from Bio import Phylo
    s = Path(path_or_str).read_text() if Path(str(path_or_str)).exists() else str(path_or_str)
    return Phylo.read(io.StringIO(s.strip()), "newick")


def midpoint(tree):
    tree.root_at_midpoint(); return tree


def parse_support(label):
    """IQ-TREE 'SH-aLRT/UFBoot' internal labels -> (sh, uf)."""
    if label is None:
        return (np.nan, np.nan)
    parts = str(label).split("/")
    try:
        if len(parts) == 2:
            return float(parts[0]), float(parts[1])
        return (np.nan, float(parts[0]))
    except ValueError:
        return (np.nan, np.nan)


def propagate_clades(tree, labels: dict, min_purity: float = 0.9) -> pd.DataFrame:
    """Assign Shah-2023 sub-clades to unlabelled tips.

    For each sub-clade, the largest clade in which >= min_purity of the labelled
    descendants belonging to it (and >=80 % of that sub-clade's labelled members) is
    taken as the sub-clade node; unlabelled tips inside the smallest such node inherit
    the label. Conflicts and unresolved tips are reported.
    """
    tips = [t.name for t in tree.get_terminals()]
    lab = {t: labels.get(t) for t in tips}
    clade_nodes = {}
    for c in sorted(set(v for v in lab.values() if v)):
        members = [t for t in tips if lab[t] == c]
        best = None
        for node in tree.find_clades(order="postorder"):
            if node.is_terminal():
                continue
            names = [x.name for x in node.get_terminals()]
            labd = [lab[n] for n in names if lab[n]]
            if not labd:
                continue
            purity = sum(1 for x in labd if x == c) / len(labd)
            recall = sum(1 for x in labd if x == c) / len(members)
            if purity >= min_purity and recall >= 0.8:
                if best is None or len(names) > best[1]:       # largest pure ancestor
                    best = (node, len(names), purity, recall)
        if best:
            clade_nodes[c] = best
    rows = []
    for t in tips:
        inferred, how = lab[t], "Shah2023"
        if not inferred:
            cands = [(c, b[1]) for c, b in clade_nodes.items() if t in {x.name for x in b[0].get_terminals()}]
            if cands:
                inferred, how = sorted(cands, key=lambda x: x[1])[0][0], "propagated"
            else:
                inferred, how = "unassigned", "none"
        rows.append(dict(tip=t, shah_label=lab[t] or "", clade=inferred, source=how))
    df = pd.DataFrame(rows)
    df.attrs["clade_nodes"] = {c: dict(size=b[1], purity=b[2], recall=b[3]) for c, b in clade_nodes.items()}
    return df


# ------------------------------------------------------------------ jplace
def parse_jplace(path) -> tuple[object, pd.DataFrame]:
    """Return (tree with edge numbers in .comment, placements DataFrame)."""
    from Bio import Phylo
    J = json.load(open(path))
    nwk = J["tree"]
    nwk = re.sub(r"\{(\d+)\}", r"[&edge=\1]", nwk)          # edge ids -> comments
    tree = Phylo.read(io.StringIO(nwk), "newick", comments_are_confidence=False)
    for node in tree.find_clades():
        m = re.search(r"edge=(\d+)", node.comment or "")
        node.edge = int(m.group(1)) if m else None
    f = J["fields"]
    rows = []
    for pl in J["placements"]:
        names = [n if isinstance(n, str) else n[0] for n in pl.get("n", pl.get("nm", []))]
        for p in pl["p"]:
            d = dict(zip(f, p))
            for nm in names:
                rows.append(dict(query=nm, **d))
    return tree, pd.DataFrame(rows)


def edge_to_tips(tree) -> dict:
    return {n.edge: [t.name for t in n.get_terminals()] for n in tree.find_clades() if getattr(n, "edge", None) is not None}


def placement_clade_summary(tree, pl: pd.DataFrame, tip_clade: dict) -> pd.DataFrame:
    e2t = edge_to_tips(tree)
    def edge_clade(e):
        cl = [tip_clade.get(t, "unassigned") for t in e2t.get(e, [])]
        cl = [c for c in cl if c != "unassigned"] or ["unassigned"]
        cnt = collections.Counter(cl)
        top, k = cnt.most_common(1)[0]
        return top if k / sum(cnt.values()) >= 0.5 else "/".join(sorted(cnt))
    pl = pl.copy()
    pl["clade"] = pl.edge_num.map(edge_clade)
    pl["n_tips_below"] = pl.edge_num.map(lambda e: len(e2t.get(e, [])))
    summ = (pl.groupby(["query", "clade"]).like_weight_ratio.sum().reset_index()
              .sort_values(["query", "like_weight_ratio"], ascending=[True, False]))
    best = summ.groupby("query").head(1).rename(columns={"clade": "assigned_clade", "like_weight_ratio": "clade_LWR"})
    top = pl.sort_values("like_weight_ratio", ascending=False).groupby("query").head(1)
    best = best.merge(top[["query", "edge_num", "like_weight_ratio", "n_tips_below"]].rename(
        columns={"like_weight_ratio": "best_edge_LWR"}), on="query")
    best["nearest_tips_below_best_edge"] = best.edge_num.map(lambda e: ";".join(e2t.get(e, [])[:6]))
    return pl, summ, best


# ------------------------------------------------------------------ haplotypes / diversity
def haplotypes(aln: dict, names: list) -> tuple[dict, pd.DataFrame]:
    """Collapse sequences to haplotypes over columns with unambiguous bases in all names."""
    A = np.array([list(aln[n].upper()) for n in names])
    ok = np.isin(A, list("ACGT")).all(0)
    B = A[:, ok]
    keys = ["".join(r) for r in B]
    hap_of, haps = {}, collections.OrderedDict()
    for n, k in zip(names, keys):
        if k not in haps:
            haps[k] = f"H{len(haps) + 1:02d}"
        hap_of[n] = haps[k]
    n = len(names); counts = collections.Counter(hap_of.values())
    hd = n / (n - 1) * (1 - sum((c / n) ** 2 for c in counts.values())) if n > 1 else float("nan")
    S = int(sum(len(set(B[:, j])) > 1 for j in range(B.shape[1])))
    pair = [np.mean(B[i] != B[j]) for i, j in itertools.combinations(range(n), 2)]
    pi = float(np.mean(pair)) if pair else float("nan")
    seq_by_hap = {h: k for k, h in haps.items()}
    stats = pd.DataFrame([dict(n_seq=n, n_sites_compared=int(ok.sum()), segregating_sites=S,
                               n_haplotypes=len(haps), Hd=hd, pi=pi)])
    return dict(hap_of=hap_of, seq_by_hap=seq_by_hap, counts=counts), stats


def haplotype_msn(seq_by_hap: dict):
    """Minimum spanning network (union of all MSTs, Bandelt et al. 1999 style) via networkx."""
    import networkx as nx
    hs = list(seq_by_hap)
    d = {(a, b): sum(x != y for x, y in zip(seq_by_hap[a], seq_by_hap[b])) for a, b in itertools.combinations(hs, 2)}
    G = nx.Graph(); G.add_nodes_from(hs)
    # union of MSTs: add all edges of weight w that connect different components of graph with edges < w
    for w in sorted(set(d.values())):
        H = nx.Graph(); H.add_nodes_from(hs); H.add_edges_from([e for e, v in d.items() if v < w])
        comp = {n: i for i, c in enumerate(nx.connected_components(H)) for n in c}
        for (a, b), v in d.items():
            if v == w and comp[a] != comp[b]:
                G.add_edge(a, b, weight=v)
    return G


# ------------------------------------------------------------------ phylogeny-trait association
class _TreeIndex:
    """Pre-computed postorder structure for fast repeated trait statistics."""
    def __init__(self, tree, tips_with_state):
        keep = set(tips_with_state)
        self.tip_names = [t.name for t in tree.get_terminals() if t.name in keep]
        pos = {n: i for i, n in enumerate(self.tip_names)}
        self.nodes = []          # (children_ids, tip_index or -1)
        self.desc = []           # tip indices below each node
        def rec(n):
            if n.is_terminal():
                if n.name not in pos:
                    return None
                nid = len(self.nodes); self.nodes.append(([], pos[n.name])); self.desc.append(np.array([pos[n.name]]))
                return nid
            ch = [c for c in (rec(c) for c in n.clades) if c is not None]
            if not ch:
                return None
            if len(ch) == 1:
                return ch[0]
            nid = len(self.nodes); self.nodes.append((ch, -1))
            self.desc.append(np.concatenate([self.desc[c] for c in ch]))
            return nid
        import sys
        sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
        rec(tree.root)
        self.internal = [k for k, (ch, t) in enumerate(self.nodes) if t < 0]

    def fitch(self, codes):
        sets = [None] * len(self.nodes); score = 0
        for k, (ch, t) in enumerate(self.nodes):     # nodes are stored in postorder
            if t >= 0:
                sets[k] = {codes[t]}
            else:
                cnt = collections.Counter(itertools.chain.from_iterable(sets[c] for c in ch))
                m = max(cnt.values()); score += len(ch) - m
                sets[k] = {x for x, v in cnt.items() if v == m}
        return score

    def ai_mc(self, codes, k_states):
        ai = 0.0; mc = np.zeros(k_states, dtype=int)
        for nd in range(len(self.nodes)):
            d = self.desc[nd]; cnt = np.bincount(codes[d], minlength=k_states)
            n = len(d); f = cnt.max() / n
            if n >= 2 and self.nodes[nd][1] < 0:
                ai += (1 - f) / (2.0 ** (n - 1)) if n < 1000 else 0.0
            if f == 1.0:
                s = cnt.argmax(); mc[s] = max(mc[s], n)
        return ai, mc


def trait_association(tree, state: dict, n_perm: int = 1000, seed: int = SEED, mc_values=None) -> pd.DataFrame:
    """BaTS-like tests (Parker et al. 2008) on a single (ML) tree using tip-label permutation:
    parsimony score (PS), association index (AI) and maximum monophyletic clade (MC)."""
    rng = np.random.default_rng(seed)
    T = _TreeIndex(tree, [k for k, v in state.items() if v is not None])
    cats = sorted(set(state[t] for t in T.tip_names)); ci = {c: i for i, c in enumerate(cats)}
    codes = np.array([ci[state[t]] for t in T.tip_names])
    mc_values = mc_values or [c for c, _ in collections.Counter(state[t] for t in T.tip_names).most_common(5)]
    def stats(cd):
        ai, mc = T.ai_mc(cd, len(cats)); r = dict(PS=T.fitch(cd), AI=ai)
        for v in mc_values:
            r[f"MC[{v}]"] = int(mc[ci[v]])
        return r
    obs = stats(codes); null = collections.defaultdict(list)
    for _ in range(n_perm):
        for k, v in stats(rng.permutation(codes)).items():
            null[k].append(v)
    rows = []
    for k, o in obs.items():
        nl = np.array(null[k])
        p = (np.sum(nl <= o) + 1) / (n_perm + 1) if k in ("PS", "AI") else (np.sum(nl >= o) + 1) / (n_perm + 1)
        rows.append(dict(statistic=k, observed=o, null_mean=nl.mean(), null_2_5=np.percentile(nl, 2.5),
                         null_97_5=np.percentile(nl, 97.5), p=p, n_tips=len(codes), n_states=len(cats)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ Rep protein
REP_MOTIFS = {  # searched in the translated amplicon; positions reported relative to the amplicon
    "RCR motif II (HLQG)": r"HLQG",
    "RCR motif III (YxxK)": r"Y.{2}K",
    "SF3 Walker A (GxxGxGKS/T)": r"G..G.GK[ST]",
    "SF3 Walker B (hhhDD)": r"[ILVFM]{3}DD",
    "SF3 motif C (hTSN)": r"[ILVFM]T[ST]N",
}


def translate(s: str) -> str:
    from Bio.Seq import Seq
    s = s.replace("-", "")
    return str(Seq(s[: len(s) // 3 * 3]).translate())


def motif_table(seqs: dict) -> pd.DataFrame:
    rows = []
    for n, s in seqs.items():
        aa = translate(s)
        r = dict(sequence=n, aa_length=len(aa), stop_codons=aa.count("*"))
        for k, p in REP_MOTIFS.items():
            m = re.search(p, aa)
            r[k] = f"{m.group()}@{m.start() + 1}" if m else "absent"
        rows.append(r)
    return pd.DataFrame(rows)


def codon_alignment(nt: dict, mafft_bin="mafft", workdir: Path = Path("."), log=None) -> dict:
    """Translate (frame 0), align proteins with MAFFT L-INS-i, thread nucleotides back (pal2nal-like)."""
    prot = collections.OrderedDict((k, translate(v).replace("*", "X")) for k, v in nt.items())
    write_fasta(prot, workdir / "codon_input.faa")
    run(f"{mafft_bin} --localpair --maxiterate 1000 --quiet {workdir/'codon_input.faa'} > {workdir/'codon_aln.faa'}", log)
    paln = read_fasta(workdir / "codon_aln.faa")
    out = collections.OrderedDict()
    for k, pa in paln.items():
        s = nt[k].replace("-", ""); i = 0; codons = []
        for c in pa:
            if c == "-":
                codons.append("---")
            else:
                codons.append(s[i:i + 3]); i += 3
        out[k] = "".join(codons)
    return out


def hyphy_json_sites(path, method: str) -> pd.DataFrame:
    J = json.load(open(path))
    mle = J["MLE"]; hdr = [h[0] for h in mle["headers"]]
    content = mle["content"]["0"] if isinstance(mle["content"], dict) else mle["content"]
    df = pd.DataFrame(content, columns=hdr)
    df.insert(0, "codon", range(1, len(df) + 1)); df.insert(0, "method", method)
    return df
