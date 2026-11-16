"""Publication figures. Every host / species label is resolved through pbfd_names,
so re-running these functions after filling FINAL_* columns updates all figures.
"""
from __future__ import annotations
from pathlib import Path
import collections, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from pbfd_names import species_label, host_label, kerman_tip_label, load_species_groups, load_host_mapping, DAGGER

# Okabe-Ito colour-blind-safe palette
OI = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
CLADE_COL = {"GI-a": "#0072B2", "GI-b": "#E69F00", "GI-c": "#009E73", "GI-d": "#CC79A7", "GI-e": "#56B4E9",
             "GI-f": "#D55E00", "GII-a": "#F0E442", "GII-b": "#7F7F7F", "unassigned": "#DDDDDD"}
BAND_COL = {"None": "#DDDDDD", "Low": "#9ecae1", "Moderate": "#3182bd", "High": "#08306b"}
KERMAN_COL, IRAN_COL = "#D55E00", "#CC79A7"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 600, "pdf.fonttype": 42, "svg.fonttype": "none"})


def save(fig, out: Path, name: str, exts=("pdf", "png", "tiff"), dpi=600):
    out.mkdir(parents=True, exist_ok=True)
    for ext in exts:
        kw = dict(pil_kwargs={"compression": "tiff_lzw"}) if ext == "tiff" else {}
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight", dpi=dpi, **kw)
    plt.close(fig)


def pending_note(ax_or_fig):
    return  # final names applied; no pending-name footnote
    txt = f"{DAGGER} name pending verification (placeholder / suggested)"
    (ax_or_fig.text(0.0, -0.02, txt, fontsize=6, color="0.4", transform=ax_or_fig.transFigure, va="top")
     if hasattr(ax_or_fig, "transFigure") else None)


# ------------------------------------------------------------------ Fig 1: prevalence
def fig_prevalence(prev: dict, meta_dir: Path, out: Path):
    groups = load_species_groups(meta_dir)
    sp = prev["by_species"].copy()
    sp["label"] = sp.group.map(lambda g: species_label(g, groups))
    sp = sp.sort_values("prevalence")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[3, 1.1], wspace=0.9))
    y = np.arange(len(sp))
    a.errorbar(sp.prevalence * 100, y, xerr=[(sp.prevalence - sp.ci_low) * 100, (sp.ci_high - sp.prevalence) * 100],
               fmt="o", color=OI[0], ecolor="0.5", capsize=2, ms=4)
    ov = prev["overall"].iloc[0]
    a.axvline(ov.prevalence * 100, ls="--", color="0.6", lw=0.8)
    a.set_yticks(y); a.set_yticklabels([f"{l}  (n={n})" for l, n in zip(sp.label, sp.n)], fontsize=7)
    a.set_xlabel("PCR-positive birds, % (Wilson 95% CI)"); a.set_xlim(0, 100)
    a.set_title("A  Prevalence by host group", loc="left", fontweight="bold")
    sx = prev["by_sex"]
    b.bar(sx.group, sx.prevalence * 100, color=[OI[3], OI[4]][: len(sx)],
          yerr=[(sx.prevalence - sx.ci_low) * 100, (sx.ci_high - sx.prevalence) * 100], capsize=3)
    for i, r in enumerate(sx.itertuples()):
        b.text(i, 2, f"{r.positive}/{r.n}", ha="center", color="white", fontsize=7)
    b.set_ylim(0, 100); b.set_ylabel("PCR-positive, %"); b.set_title("B  By sex", loc="left", fontweight="bold")
    pending_note(fig); save(fig, out, "Fig1_prevalence")


# ------------------------------------------------------------------ Fig 2: alluvial
def fig_alluvial(birds: pd.DataFrame, meta_dir: Path, out: Path):
    """Three-axis alluvial: host group -> PCR band category -> sex."""
    groups = load_species_groups(meta_dir)
    birds = birds.copy(); birds["lab"] = birds.species_group.map(lambda g: species_label(g, groups))
    axes = ["lab", "band_label", "sex"]
    order = {"lab": list(birds.lab.value_counts().index), "band_label": ["None", "Low", "Moderate", "High"],
             "sex": ["female", "male"]}
    N, gap = len(birds), 4
    pos = {}
    for ai, ax in enumerate(axes):
        y = 0
        for v in order[ax]:
            n = int((birds[ax].astype(str) == v).sum())
            pos[(ai, v)] = [y, y + n, y]  # start, end, cursor
            y += n + gap
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    band_of = birds.band_label.astype(str)
    for (i, a1), (j, a2) in [((0, "lab"), (1, "band_label")), ((1, "band_label"), (2, "sex"))]:
        flows = birds.groupby([birds[a1].astype(str), birds[a2].astype(str)], observed=True).size()
        for v1 in order[a1]:
            for v2 in order[a2]:
                n = int(flows.get((v1, v2), 0))
                if not n:
                    continue
                s = pos[(i, v1)]; t = pos[(j, v2)]
                y0, y1 = s[2], t[2]; s[2] += n; t[2] += n
                col = BAND_COL[v1 if a1 == "band_label" else v2]
                xs = np.linspace(i + 0.06, j - 0.06, 50); w = (1 - np.cos(np.linspace(0, np.pi, 50))) / 2
                top = y0 + (y1 - y0) * w; bot = top + n
                ax.fill_between(xs, top, bot, color=col, alpha=0.55, lw=0)
        for v in order[a1] if i == 0 else []:
            pass
    for ai, a_ in enumerate(axes):
        for v in order[a_]:
            s, e, _ = pos[(ai, v)]
            if e == s:
                continue
            ax.add_patch(Rectangle((ai - 0.06, s), 0.12, e - s, color="0.25" if a_ != "band_label" else BAND_COL[v], ec="k", lw=0.3))
            lab = f"{v} ({e - s})"
            ha = "right" if ai == 0 else "left"; x = ai - 0.08 if ai == 0 else ai + 0.08
            ax.text(x, (s + e) / 2, lab, ha=ha, va="center", fontsize=6.5)
    ax.set_xlim(-1.3, 2.6); ax.invert_yaxis(); ax.axis("off")
    for ai, t in enumerate(["Host group", "PCR band intensity", "Sex"]):
        ax.text(ai, -6, t, ha="center", fontweight="bold")
    pending_note(fig); save(fig, out, "Fig2_alluvial_host_band_sex")


# ------------------------------------------------------------------ tree drawing
def nice_scale(v):
    if v <= 0:
        return 0.01
    e = 10 ** math.floor(math.log10(v))
    return max(m * e for m in (1, 2, 5) if m * e <= v)


def _layout(tree):
    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}
    x = {}
    def rec(n, d):
        x[n] = d
        for c in n.clades:
            rec(c, d + (c.branch_length or 0))
    rec(tree.root, 0)
    for n in tree.get_nonterminals(order="postorder"):
        y[n] = (y[n.clades[0]] + y[n.clades[-1]]) / 2
    return x, y


def draw_tree(ax, tree, tip_label=lambda n: n, tip_colour=lambda n: "k", clade_of=None, support=True,
              fontsize=5, lw=0.5, strip_x=None):
    from pbfd_phylo import parse_support
    x, y = _layout(tree)
    xmax = max(x.values())
    for n in tree.find_clades():
        for c in n.clades:
            ax.plot([x[n], x[n]], [y[n], y[c]], color="k", lw=lw)
            ax.plot([x[n], x[c]], [y[c], y[c]], color="k", lw=lw)
        if support and not n.is_terminal() and n.clades:
            sh, uf = parse_support(n.name if n.name else n.confidence)
            if not np.isnan(uf) and uf >= 95 and (np.isnan(sh) or sh >= 80) and len(n.get_terminals()) > 2:
                ax.plot(x[n], y[n], "o", ms=2.2, color="k")
    for t in tree.get_terminals():
        ax.text(x[t] + xmax * 0.01, y[t], tip_label(t.name), fontsize=fontsize, va="center", color=tip_colour(t.name))
    if clade_of is not None:
        sx = strip_x if strip_x is not None else xmax * 1.9
        for t in tree.get_terminals():
            ax.add_patch(Rectangle((sx, y[t] - 0.5), xmax * 0.04, 1, color=CLADE_COL.get(clade_of.get(t.name, "unassigned"), "#ddd"), lw=0))
    ax.axis("off")
    # scale bar
    sb = nice_scale(xmax / 5)
    ny = len(tree.get_terminals())
    ax.plot([0, sb], [ny + 0.3] * 2, color="k", lw=1, clip_on=False)
    ax.text(sb / 2, ny + 0.6, f"{sb:g} subst./site", ha="center", va="top", fontsize=6)
    ax.set_ylim(ny + 1.5, -1)
    return x, y


def subset_tree(tree, keep: set):
    import copy
    t = copy.deepcopy(tree)
    for tip in [x for x in t.get_terminals() if x.name not in keep]:
        t.prune(tip)
    return t


def _ref_label(name, refmd):
    r = refmd.get(name, {})
    year = r.get("year_or_submit", "")
    year = "" if pd.isna(year) or year == "" else str(int(year))
    host = r.get("host", "")
    host = "" if (host is None or (isinstance(host, float) and pd.isna(host))) else str(host)
    host = host if host else "host not recorded"
    return f"{name} {r.get('country_simple', '')} {year} | {host}".strip()


def fig_tree_focus(tree, keep, tip_clade: dict, refmd: dict, meta_dir: Path, out: Path, name="Fig3_ML_tree_focus",
                   title="Maximum-likelihood tree (ORF1/Rep amplicon)"):
    hosts = load_host_mapping(meta_dir)
    t = subset_tree(tree, set(keep))
    t.ladderize()
    n = len(t.get_terminals())
    fig, ax = plt.subplots(figsize=(7.2, max(4, n * 0.13)))
    def lab(nm):
        return kerman_tip_label(nm, hosts) if "_Kerman_" in nm else _ref_label(nm, refmd)
    def col(nm):
        if "_Kerman_" in nm:
            return KERMAN_COL
        return IRAN_COL if refmd.get(nm, {}).get("country_simple") == "Iran" else "k"
    draw_tree(ax, t, lab, col, tip_clade, fontsize=5.5)
    handles = [Rectangle((0, 0), 1, 1, color=c) for k, c in CLADE_COL.items() if k in set(tip_clade.get(x.name) for x in t.get_terminals())]
    labels = [k for k in CLADE_COL if k in set(tip_clade.get(x.name) for x in t.get_terminals())]
    ax.legend(handles, labels, title="Sub-clade (Shah et al. 2023)", loc="upper center", bbox_to_anchor=(0.5, -0.005), fontsize=6, title_fontsize=6, frameon=False, ncol=4)
    ax.set_title(title + f"\n● UFBoot ≥95 & SH-aLRT ≥80; orange = this study; pink = other Iranian strains", fontsize=7, loc="left")
    save(fig, out, name)


def fig_tree_full(tree, tip_clade, refmd, meta_dir, out, name="FigS1_ML_tree_full"):
    hosts = load_host_mapping(meta_dir)
    tree.ladderize()
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(8.3, n * 0.075 + 1))
    lab = lambda nm: kerman_tip_label(nm, hosts) if "_Kerman_" in nm else _ref_label(nm, refmd)
    col = lambda nm: KERMAN_COL if "_Kerman_" in nm else (IRAN_COL if refmd.get(nm, {}).get("country_simple") == "Iran" else "k")
    draw_tree(ax, tree, lab, col, tip_clade, fontsize=3, lw=0.3)
    save(fig, out, name, exts=("pdf", "png"), dpi=200)


# ------------------------------------------------------------------ Fig 4: placement
def fig_placement(summ: pd.DataFrame, meta_dir: Path, out: Path):
    hosts = load_host_mapping(meta_dir)
    piv = summ.pivot_table(index="query", columns="clade", values="like_weight_ratio", aggfunc="sum").fillna(0)
    piv = piv.loc[sorted(piv.index)]
    fig, ax = plt.subplots(figsize=(7.2, 0.35 * len(piv) + 1.2))
    left = np.zeros(len(piv))
    for c in piv.columns:
        ax.barh(range(len(piv)), piv[c], left=left, color=CLADE_COL.get(c, "#999"), label=c, ec="white", lw=0.3)
        left += piv[c].values
    ax.set_yticks(range(len(piv))); ax.set_yticklabels([kerman_tip_label(q, hosts) for q in piv.index], fontsize=6.5)
    ax.set_xlim(0, 1); ax.set_xlabel("Cumulative likelihood weight ratio (EPA-ng) on full-genome backbone")
    ax.invert_yaxis(); ax.legend(ncol=4, fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    save(fig, out, "Fig4_placement_backbone")


# ------------------------------------------------------------------ Fig 5: identity + network
def fig_identity_network(D: pd.DataFrame, G, hap_of: dict, hap_counts, node_attr: dict, meta_dir: Path, out: Path):
    import networkx as nx
    hosts = load_host_mapping(meta_dir)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.8), gridspec_kw=dict(width_ratios=[1.15, 1]))
    ident = (1 - D) * 100
    short = [n.split("_Kerman_")[-1] + "*" if "_Kerman_" in n else n for n in ident.index]
    im = a.imshow(ident.values, cmap="viridis", vmin=max(85, np.nanmin(ident.values)), vmax=100)
    a.set_xticks(range(len(short))); a.set_xticklabels(short, rotation=90, fontsize=5)
    a.set_yticks(range(len(short))); a.set_yticklabels(short, fontsize=5)
    for i in range(len(short)):
        for j in range(len(short)):
            if len(short) <= 22:
                a.text(j, i, f"{ident.values[i, j]:.1f}", ha="center", va="center", fontsize=3.4,
                       color="w" if ident.values[i, j] < 97 else "k")
    fig.colorbar(im, ax=a, fraction=0.04, label="nt identity, %")
    a.set_title("A  Pairwise identity (* = this study)", loc="left", fontweight="bold")
    pos = nx.kamada_kawai_layout(G, weight="weight") if len(G) > 1 else {n: (0, 0) for n in G}
    for u, v, d in G.edges(data=True):
        b.plot(*zip(pos[u], pos[v]), color="0.5", lw=0.7, zorder=1)
        mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
        b.text(mx, my, str(d["weight"]), fontsize=5, color="0.3", ha="center", va="center",
               bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none"))
    for h in G.nodes:
        members = [n for n, hh in hap_of.items() if hh == h]
        groups = collections.Counter(node_attr.get(m, "other") for m in members)
        r = 0.035 + 0.02 * math.sqrt(len(members))
        start = 0
        cols = {"this study (Kerman)": KERMAN_COL, "Iran (other)": IRAN_COL}
        for gname, k in groups.items():
            ang = 360 * k / len(members)
            from matplotlib.patches import Wedge
            b.add_patch(Wedge(pos[h], r, start, start + ang, color=cols.get(gname, "0.6"), zorder=2, ec="k", lw=0.3))
            start += ang
        b.text(pos[h][0], pos[h][1] - r - 0.04, h, fontsize=5, ha="center")
    b.set_aspect("equal"); b.axis("off")
    b.set_title("B  Minimum-spanning haplotype network\n(numbers = substitutions)", loc="left", fontweight="bold", fontsize=7)
    from matplotlib.patches import Patch
    b.legend([Patch(color=KERMAN_COL), Patch(color=IRAN_COL), Patch(color="0.6")],
             ["this study (Kerman)", "Iran (other)", "nearest non-Iranian"], fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=3)
    save(fig, out, "Fig5_identity_haplotype_network")


# ------------------------------------------------------------------ supplementary
def fig_sanger_qc(qc: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(5, 2.6))
    ax.bar(qc.voucher, qc.q20_fraction * 100, color=OI[0])
    for i, r in enumerate(qc.itertuples()):
        ax.text(i, r.q20_fraction * 100 + 1, f"{r.mean_phred:.0f}", ha="center", fontsize=6)
    ax.set_ylabel("Bases with Phred ≥20, %"); ax.set_ylim(0, 105)
    ax.set_title("Sanger read quality (label = mean Phred)", loc="left", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    save(fig, out, "FigS2_sanger_qc")


def fig_selection(sites: pd.DataFrame, motifs: dict, out: Path):
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    for meth, mk, col in [("FEL", "o", OI[0]), ("MEME", "^", OI[5])]:
        s = sites[sites.method == meth]
        pcol = next((c for c in s.columns if c.lower().startswith("p-value")), None)
        if pcol is None:
            continue
        ax.scatter(s.codon, -np.log10(s[pcol].astype(float).clip(lower=1e-6)), s=8, marker=mk, color=col, label=meth)
    ax.axhline(-np.log10(0.1), ls="--", color="0.5", lw=0.7)
    for k, (a0, a1) in motifs.items():
        ax.axvspan(a0, a1, color="0.9", zorder=0); ax.text((a0 + a1) / 2, ax.get_ylim()[1] * 0.95, k.split(" (")[0], fontsize=5, ha="center", rotation=90, va="top")
    ax.set_xlabel("Codon position in alignment"); ax.set_ylabel("-log10 p"); ax.legend(fontsize=6, frameon=False)
    save(fig, out, "FigS3_selection_sites")
