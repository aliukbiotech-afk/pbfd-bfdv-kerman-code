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

# ---------------------------------------------------------------------------
# UNIFIED COLOUR TOKENS - one meaning per colour, identical in every main and
# supplementary figure.  (Documented in FIGURE_COLOUR_KEY.md.)
#   clade   : Shah et al. 2023 sub-clades (trees, placement)
#   origin  : provenance of a sequence (this study / other Iran / non-Iranian)
#   host    : host group of a bird or isolate (survey and QC figures)
#   result  : PCR outcome;  sex : neutral slate tones (never a host/clade colour)
# ---------------------------------------------------------------------------
CLADE_COL = {"GI-a": "#0072B2", "GI-b": "#56B4E9", "GI-c": "#009E73", "GI-d": "#D9B300", "GI-e": "#332288",
             "GI-f": "#882255", "GII-a": "#8C8C8C", "GII-b": "#4D4D4D", "unassigned": "#DDDDDD"}
KERMAN_COL, IRAN_COL, OTHER_COL = "#D55E00", "#CC79A7", "#707070"
HOST_COL = {"SP01": "#264653", "SP02": "#2A9D8F", "SP03": "#8AB17D", "SP04": "#E9C46A",
            "SP05": "#F4A261", "SP06": "#E76F51", "SP07": "#9A9A9A"}
RESULT_COL = {"Positive": "#2B2B2B", "Negative": "#D0D0D0"}
SEX_COL = {"female": "#5B6770", "male": "#A9B4BC"}
OI = [HOST_COL["SP01"], KERMAN_COL, CLADE_COL["GI-c"], IRAN_COL]          # legacy aliases
INK, GRID = "#222222", "#E6E6E6"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#444444", "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": "#444444", "ytick.color": "#444444",
                     "axes.titleweight": "bold", "axes.titlesize": 8.5, "savefig.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


def panel(ax, letter, title, x=0.0, y=1.04):
    ax.text(x, y, f"{letter}", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")
    ax.text(x + 0.015, y + 0.003, title, transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="bottom", ha="left")


def host_group_of(host_name, groups):
    r = groups.loc[groups.FINAL_SPECIES == host_name]
    return r.group_id.iloc[0] if len(r) else "SP07"


def iso_label(voucher, hosts):
    """Specimen label = host species only (isolate IDs live in the tables)."""
    return host_label(voucher, hosts)


def save(fig, out: Path, name: str, exts=("pdf", "png", "tiff"), dpi=600):
    out.mkdir(parents=True, exist_ok=True)
    for ext in exts:
        kw = dict(pil_kwargs={"compression": "tiff_lzw"}) if ext == "tiff" else {}
        for attempt in range(8):                      # a sync client may lock files briefly
            try:
                fig.savefig(out / f"{name}.{ext}", bbox_inches="tight", dpi=dpi, **kw); break
            except OSError:
                if attempt == 7:
                    raise
                __import__("time").sleep(1.5)
    plt.close(fig)


# ------------------------------------------------------------------ Fig 1: prevalence
def fig_prevalence(prev: dict, meta_dir: Path, out: Path):
    groups = load_species_groups(meta_dir)
    sp = prev["by_species"].copy()
    sp["label"] = sp.group.map(lambda g: species_label(g, groups))
    sp = sp.sort_values("prevalence")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(width_ratios=[3, 1.1], wspace=0.62))
    y = np.arange(len(sp))
    a.hlines(y, sp.ci_low * 100, sp.ci_high * 100, color="0.6", lw=1.2, zorder=1)
    a.scatter(sp.prevalence * 100, y, s=34, c=[HOST_COL[g] for g in sp.group], ec=INK, lw=0.5, zorder=2)
    ov = prev["overall"].iloc[0]
    a.axvline(ov.prevalence * 100, ls="--", color="0.55", lw=0.8)
    a.text(ov.prevalence * 100 + 1, len(sp) - 0.45, f"all birds {ov.prevalence*100:.1f}%", fontsize=6.5, color="0.35", va="bottom")
    a.set_yticks(y); a.set_yticklabels([f"{l}  (n={n})" for l, n in zip(sp.label, sp.n)], fontsize=7)
    a.set_xlabel("PCR-positive birds, % (Wilson 95% CI)"); a.set_xlim(0, 100); a.set_ylim(-0.6, len(sp) + 0.3)
    a.grid(axis="x", color=GRID, lw=0.6, zorder=0); a.set_axisbelow(True)
    panel(a, "A", "Prevalence by host group", x=-0.02)
    sx = prev["by_sex"]
    b.bar(sx.group, sx.prevalence * 100, color=[SEX_COL.get(g, "0.6") for g in sx.group], ec=INK, lw=0.5,
          yerr=[(sx.prevalence - sx.ci_low) * 100, (sx.ci_high - sx.prevalence) * 100], capsize=3, error_kw=dict(lw=0.8, ecolor="0.3"))
    for i, r in enumerate(sx.itertuples()):
        b.text(i, 3, f"{r.positive}/{r.n}", ha="center", color="white" if r.group == "female" else INK, fontsize=7)
    b.set_ylim(0, 100); b.set_ylabel("PCR-positive, %"); panel(b, "B", "By sex", x=-0.02)
    save(fig, out, "Fig1_prevalence")


# ------------------------------------------------------------------ Fig 2: alluvial (host -> PCR result -> sex)
def fig_alluvial(birds: pd.DataFrame, meta_dir: Path, out: Path, note: str | None = None):
    """Purpose: show, for all 136 birds, how PCR-positive birds are distributed across host groups and sexes.
    Every column sums to the same 136 birds and is drawn at the same total height.
    PCR band intensity is deliberately not shown (semi-quantitative; reported in the text only)."""
    groups = load_species_groups(meta_dir)
    d = birds.copy()
    d["host"] = d.species_group
    d["result"] = np.where(d.positive == 1, "Positive", "Negative")
    cols = ["host", "result", "sex"]
    hosts_order = list(d.host.value_counts().index)
    order = {"host": hosts_order, "result": ["Positive", "Negative"], "sex": ["female", "male"]}
    N, SPAN = len(d), 136 + 48
    pos = {}
    for ci, c in enumerate(cols):
        k = len(order[c]); gap = (SPAN - N) / (k - 1); y = 0.0
        for v in order[c]:
            n = int((d[c] == v).sum()); pos[(ci, v)] = [y, y + n, y, y]; y += n + gap
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    w = 0.07

    def flow(ci, v1, v2, n, col):
        s = pos[(ci, v1)]; t = pos[(ci + 1, v2)]
        y0, y1 = s[2], t[3]; s[2] += n; t[3] += n
        xs = np.linspace(ci + w, ci + 1 - w, 60); k = (1 - np.cos(np.linspace(0, np.pi, 60))) / 2
        top = y0 + (y1 - y0) * k
        ax.fill_between(xs, top, top + n, color=col, alpha=0.6, lw=0)

    for r in order["result"]:                      # positives first so they sit at the top of each host block
        pass
    for h in hosts_order:
        for r in order["result"]:
            n = int(((d.host == h) & (d.result == r)).sum())
            if n:
                flow(0, h, r, n, HOST_COL[h])
    for r in order["result"]:
        for sx in order["sex"]:
            n = int(((d.result == r) & (d.sex == sx)).sum())
            if n:
                flow(1, r, sx, n, RESULT_COL[r])
    # host labels: one line each, spread so none touch, with thin leader lines
    ys = {v: (pos[(0, v)][0] + pos[(0, v)][1]) / 2 for v in hosts_order}
    lab_y = dict(ys); min_sp = 8.5
    for _ in range(200):
        moved = False
        for v1, v2 in zip(hosts_order[:-1], hosts_order[1:]):
            if lab_y[v2] - lab_y[v1] < min_sp:
                m = (lab_y[v1] + lab_y[v2]) / 2
                lab_y[v1] = m - min_sp / 2; lab_y[v2] = m + min_sp / 2; moved = True
        if not moved:
            break
    for ci, c in enumerate(cols):
        for v in order[c]:
            s, e = pos[(ci, v)][:2]
            colr = HOST_COL[v] if c == "host" else RESULT_COL[v] if c == "result" else SEX_COL[v]
            ax.add_patch(Rectangle((ci - w, s), 2 * w, e - s, color=colr, ec=INK, lw=0.4, zorder=3))
            n_v = int(e - s)
            if c == "host":
                npos = int(((d.host == v) & (d.result == "Positive")).sum())
                ax.plot([ci - w, ci - w - 0.05, ci - w - 0.09], [ys[v], lab_y[v], lab_y[v]], color="0.55", lw=0.5)
                ax.text(ci - w - 0.11, lab_y[v], f"{species_label(v, groups)}   {npos}/{n_v} positive", ha="right", va="center", fontsize=6.6)
            else:
                lab = f"{v.capitalize()}\n{n_v} birds"
                ax.text(ci + w + 0.03, (s + e) / 2, lab, ha="left", va="center", fontsize=7, linespacing=1.15,
                        bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.2), zorder=4)
    ax.set_xlim(-2.25, 2.55); ax.set_ylim(SPAN + 8, -14); ax.axis("off")
    for ci, t in enumerate(["Host group", "PCR result", "Sex"]):
        ax.text(ci, -9, t, ha="center", fontweight="bold", fontsize=8.5)
    if note:
        ax.text(0.5, SPAN + 5, note, ha="center", va="top", fontsize=7, color="0.25", style="italic")
    save(fig, out, "Fig2_alluvial_host_result_sex")


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


VERNACULAR_TO_SCI = {"budgerigar": "Melopsittacus undulatus", "cockatiel": "Nymphicus hollandicus",
                     "love bird": "Agapornis spp.", "jenday conure": "Aratinga jandaya",
                     "parrot": "parrot (species not specified)", "psittacine bird": "psittacine (species not specified)"}


def _ref_label(name, refmd):
    r = refmd.get(name, {})
    year = r.get("year_or_submit", "")
    year = "" if pd.isna(year) or year == "" else str(int(year))
    host = r.get("host", "")
    host = "" if (host is None or (isinstance(host, float) and pd.isna(host))) else str(host)
    host = VERNACULAR_TO_SCI.get(host.strip().lower(), host) if host else "host not recorded"
    return f"{host} · {r.get('country_simple', '')} {year}".strip()


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
    """Question answered: to which Shah et al. (2023) sub-clade does each short Kerman Rep sequence belong,
    and how sure is the phylogenetic placement?  Bar = share of the likelihood weight (LWR) falling in each clade."""
    hosts = load_host_mapping(meta_dir)
    piv = summ.pivot_table(index="query", columns="clade", values="like_weight_ratio", aggfunc="sum").fillna(0)
    piv["other / unlabelled"] = (1 - piv.sum(axis=1)).clip(lower=0)
    top = piv.drop(columns="other / unlabelled").idxmax(axis=1)
    key = pd.DataFrame({"clade": top, "lwr": piv.max(axis=1)})
    order = key.sort_values(["clade", "lwr"], ascending=[True, False]).index
    piv = piv.loc[order]
    fig, ax = plt.subplots(figsize=(7.2, 0.36 * len(piv) + 1.9))
    left = np.zeros(len(piv))
    for c in [c for c in piv.columns if c != "other / unlabelled"] + ["other / unlabelled"]:
        v = piv[c].values
        if c == "other / unlabelled":
            ax.barh(range(len(piv)), v, left=left, color="white", hatch="////", ec="#AAAAAA", lw=0.4, label=c)
        else:
            ax.barh(range(len(piv)), v, left=left, color=CLADE_COL.get(c, "#999"), label=c, ec="white", lw=0.4)
        left += v
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([iso_label(q.split("_Kerman_")[-1], hosts) for q in piv.index], fontsize=6.8, style="italic")
    for t in ax.get_yticklabels():
        pass
    ax.set_xlim(0, 1); ax.set_xlabel("Likelihood weight ratio (LWR): probability that the sequence belongs to each sub-clade", fontsize=7.5)
    ax.invert_yaxis(); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, l, ncol=len(l), fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.43, -0.15), title="Sub-clade (Shah et al. 2023)", title_fontsize=6.5)
    ax.set_title("Genotype assignment of the 10 Kerman Rep sequences by phylogenetic placement on 429 full genomes", loc="left", fontsize=8, pad=6)
    save(fig, out, "Fig4_placement_backbone")


# ------------------------------------------------------------------ Fig 5: identity + network
def _seg_cross(p1, p2, p3, p4):
    def o(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return o(p1, p2, p3) * o(p1, p2, p4) < 0 and o(p3, p4, p1) * o(p3, p4, p2) < 0


def _pt_seg(p, a, b):
    ab = b - a; t = np.clip(np.dot(p - a, ab) / (np.dot(ab, ab) + 1e-12), 0, 1)
    return np.linalg.norm(p - (a + t * ab))


def network_layout(G, radius: dict, seed=11, tries=2500, width=2.15):
    """Distance-faithful haplotype-network layout that fills a landscape canvas.
    Edge length grows with sqrt(substitutions).  Random restarts are scored on edge crossings,
    edges running through foreign nodes, node crowding and stress; the winner is relaxed so that
    no two nodes (plus label room) touch."""
    import networkx as nx, itertools
    rng = np.random.default_rng(seed)
    nodes = list(G.nodes); edges = list(G.edges)
    H = nx.Graph(); H.add_weighted_edges_from([(u, v, 0.45 + 0.38 * math.sqrt(d["weight"])) for u, v, d in G.edges(data=True)])
    sp = dict(nx.all_pairs_dijkstra_path_length(H))
    dist = {u: {v: sp[u].get(v, 4.0) for v in nodes} for u in nodes}
    best, best_score = None, 1e18
    for k in range(tries):
        init = {n: rng.normal(size=2) * 1.5 for n in nodes}
        pos = nx.kamada_kawai_layout(G, dist=dist, pos=init)
        P = {n: np.array(pos[n]) for n in nodes}
        xs = np.array([p[0] for p in P.values()]); ys = np.array([p[1] for p in P.values()])
        Q = {n: np.array([(p[0] - xs.min()) / (np.ptp(xs) or 1) * width, (p[1] - ys.min()) / (np.ptp(ys) or 1)]) for n, p in P.items()}
        cross = sum(_seg_cross(Q[a], Q[b], Q[c], Q[d_]) for (a, b), (c, d_) in itertools.combinations(edges, 2) if len({a, b, c, d_}) == 4)
        thru = sum(_pt_seg(Q[n], Q[a], Q[b]) < radius[n] + 0.07 for (a, b) in edges for n in nodes if n not in (a, b))
        close = sum(np.linalg.norm(Q[a] - Q[b]) < 0.17 for a, b in itertools.combinations(nodes, 2))
        elen = sum(np.linalg.norm(Q[a] - Q[b]) < 0.2 for a, b in edges)
        score = 100 * cross + 60 * thru + 25 * close + 10 * elen
        if score < best_score:
            best, best_score = Q, score
    P = best
    for _ in range(800):
        moved = False
        for a, b in itertools.combinations(nodes, 2):
            dv = P[a] - P[b]; dd = np.linalg.norm(dv) + 1e-9
            need = radius[a] + radius[b] + 0.21
            if dd < need:
                push = dv / dd * (need - dd) / 2
                P[a] = P[a] + push; P[b] = P[b] - push; moved = True
        if not moved:
            break
    return P


def _bbox(ax, txt):
    r = ax.figure.canvas.get_renderer()
    bb = txt.get_window_extent(r).transformed(ax.transData.inverted())
    return (bb.x0, bb.x1, bb.y0, bb.y1)


def _overlap(a, b, pad=0.004):
    return not (a[1] + pad < b[0] or a[0] - pad > b[1] or a[3] + pad < b[2] or a[2] - pad > b[3])


VERNACULAR_ABBR = {"parrot (species not specified)": "unspecified parrot", "psittacine (species not specified)": "unspecified psittacine"}


def abbr_species(name: str) -> str:
    """'Nymphicus hollandicus' -> 'N. hollandicus' (used only where space is tight: Fig. 5)."""
    name = VERNACULAR_ABBR.get(name, name)
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) == 2 and parts[1] != "spp." else name


def _host_of_tip(tip, hosts, refmd):
    if "_Kerman_" in tip:
        return host_label(tip.split("_Kerman_")[-1], hosts), "Kerman"
    r = refmd.get(tip, {}) if refmd else {}
    h = r.get("host", "")
    h = "" if (h is None or (isinstance(h, float) and pd.isna(h))) else str(h)
    h = VERNACULAR_TO_SCI.get(h.strip().lower(), h) if h else "host not recorded"
    return h, r.get("country_simple", "")


def fig_identity_network(D: pd.DataFrame, G, hap_of: dict, hap_counts, node_attr: dict, meta_dir: Path, out: Path, refmd: dict | None = None):
    """A: pairwise identity heat-map.  B: minimum-spanning haplotype network, drawn full-width below A.
    Every label is a host species (genus abbreviated for space); isolate / accession IDs are in Supplementary Tables S6 and S7."""
    from matplotlib.patches import Wedge, Patch
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    hosts = load_host_mapping(meta_dir)
    fig = plt.figure(figsize=(7.4, 9.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.86, 1.0], hspace=0.27)
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1])
    # ---- A: heat-map
    ident = (1 - D) * 100
    names = []
    for n in ident.index:
        h, c = _host_of_tip(n, hosts, refmd)
        names.append(f"{abbr_species(h)}{' *' if c == 'Kerman' else ''}")
    im = a.imshow(ident.values, cmap="viridis", vmin=max(85, np.nanmin(ident.values)), vmax=100, aspect="auto")
    a.set_xticks(range(len(names))); a.set_xticklabels(names, rotation=90, fontsize=5.4)
    a.set_yticks(range(len(names))); a.set_yticklabels(names, fontsize=5.4)
    for lab in list(a.get_xticklabels()) + list(a.get_yticklabels()):
        if "*" in lab.get_text():
            lab.set_color(KERMAN_COL); lab.set_fontweight("bold")
    cax = make_axes_locatable(a).append_axes("right", size="2.2%", pad=0.08)
    cb = fig.colorbar(im, cax=cax); cb.set_label("nucleotide identity, %", fontsize=7); cb.ax.tick_params(labelsize=6)
    a.set_title("Pairwise nucleotide identity (* = this study)", loc="left", fontsize=8.5, fontweight="bold", pad=8)
    a.text(-0.01, 1.045, "A", transform=a.transAxes, fontsize=12, fontweight="bold", ha="right", va="bottom")
    # ---- B: haplotype network
    members = {h: [n for n, hh in hap_of.items() if hh == h] for h in G.nodes}
    radius = {h: 0.034 + 0.016 * math.sqrt(len(m)) for h, m in members.items()}
    pos = network_layout(G, radius, tries=700, width=2.5)
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    b.set_xlim(min(xs) - 0.6, max(xs) + 0.6); b.set_ylim(min(ys) - 0.34, max(ys) + 0.3)
    b.set_aspect("equal", adjustable="box"); b.axis("off")
    for u, v, d in G.edges(data=True):
        b.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="0.5", lw=1.1, zorder=1, solid_capstyle="round")
    cols = {"this study (Kerman)": KERMAN_COL, "Iran (other)": IRAN_COL}
    for h in G.nodes:
        grp = collections.Counter(node_attr.get(m, "other") for m in members[h]); start = 0
        for gname, k in grp.items():
            ang = 360 * k / len(members[h])
            b.add_patch(Wedge(pos[h], radius[h], start, start + ang, color=cols.get(gname, OTHER_COL), zorder=3, ec="k", lw=0.5))
            start += ang
    fig.canvas.draw()
    node_boxes = {h: (pos[h][0] - radius[h], pos[h][0] + radius[h], pos[h][1] - radius[h], pos[h][1] + radius[h]) for h in G.nodes}
    placed = []
    for u, v, d in sorted(G.edges(data=True), key=lambda e: -e[2]["weight"]):
        tx = b.text(0, 0, str(d["weight"]), fontsize=6.2, color="0.15", ha="center", va="center", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="0.7", lw=0.3))
        bestc = None
        for t in (0.5, 0.4, 0.6, 0.32, 0.68, 0.24, 0.76):
            p = pos[u] + (pos[v] - pos[u]) * t
            tx.set_position(p); bb = _bbox(b, tx)
            pen = sum(_overlap(bb, nb, 0.01) for nb in node_boxes.values()) * 5 + sum(_overlap(bb, pb) for pb in placed) * 2
            if bestc is None or pen < bestc[0]:
                bestc = (pen, p)
            if pen == 0:
                break
        tx.set_position(bestc[1]); placed.append(_bbox(b, tx))

    def nlabel(h):
        ms = members[h]
        hh, cc = _host_of_tip(ms[0], hosts, refmd)
        base = abbr_species(hh) if "_Kerman_" in ms[0] else f"{abbr_species(hh)} ({cc})"
        return base + (f"  +{len(ms) - 1}" if len(ms) > 1 else "")

    segs = [(pos[u], pos[v]) for u, v in G.edges]
    for h in sorted(G.nodes, key=lambda x: (-len(nlabel(x)), x)):
        tx = b.text(0, 0, nlabel(h), fontsize=6.6, zorder=5, style="italic")
        bestc = None
        for k, rr in enumerate((0.014, 0.075, 0.14, 0.23)):
            r = radius[h] + rr
            for ang in range(0, 360, 20):
                ca, sa = math.cos(math.radians(ang)), math.sin(math.radians(ang))
                ha = "center" if abs(ca) < 0.35 else ("left" if ca > 0 else "right")
                va = "center" if abs(sa) < 0.35 else ("top" if sa < 0 else "bottom")
                x, y = pos[h][0] + ca * r, pos[h][1] + sa * r
                tx.set_position((x, y)); tx.set_ha(ha); tx.set_va(va); bb = _bbox(b, tx)
                pen = sum(_overlap(bb, pb, 0.015) for pb in placed) * 30
                pen += sum(_overlap(bb, nb, 0.015) for o, nb in node_boxes.items() if o != h) * 30
                pen += 4 * sum(any(_pt_seg(np.array([cx, cy]), p, q) < 0.035 for cx in np.linspace(bb[0], bb[1], 10) for cy in np.linspace(bb[2], bb[3], 3))
                               for p, q in segs)
                pen += 0.6 * k + (0.3 if ang % 90 else 0)
                if bestc is None or pen < bestc[0]:
                    bestc = (pen, x, y, ha, va, k, ca, sa)
        _, x, y, ha, va, k, ca, sa = bestc
        tx.set_position((x, y)); tx.set_ha(ha); tx.set_va(va); placed.append(_bbox(b, tx))
        if k > 0:
            b.plot([pos[h][0] + ca * radius[h], x - ca * 0.01], [pos[h][1] + sa * radius[h], y - sa * 0.01], color="0.6", lw=0.5, zorder=2)
    b.set_title("Minimum-spanning haplotype network (numbers on links = nucleotide substitutions)", loc="left", fontsize=8.5, fontweight="bold", pad=8)
    b.text(0.0, 1.04, "B", transform=b.transAxes, fontsize=12, fontweight="bold", ha="right", va="bottom")
    b.legend([Patch(color=KERMAN_COL, ec="k", lw=0.5), Patch(color=IRAN_COL, ec="k", lw=0.5), Patch(color=OTHER_COL, ec="k", lw=0.5)],
             ["this study (Kerman)", "other Iranian strains", "nearest non-Iranian strains"], fontsize=7, frameon=False,
             loc="lower center", bbox_to_anchor=(0.5, -0.03), ncol=3, handlelength=1.3)
    save(fig, out, "Fig5_identity_haplotype_network")


# ------------------------------------------------------------------ supplementary
def fig_sanger_qc(qc: pd.DataFrame, out: Path, meta_dir: Path | None = None):
    meta_dir = Path(meta_dir) if meta_dir else Path(__file__).resolve().parents[1] / "metadata"
    hosts = load_host_mapping(meta_dir); groups = load_species_groups(meta_dir)
    q = qc.copy()
    q["host"] = q.voucher.map(lambda v: host_label(v, hosts))
    q["grp"] = q.host.map(lambda h: host_group_of(h, groups))
    q = q.sort_values("q20_fraction", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    y = np.arange(len(q))
    ax.barh(y, q.q20_fraction * 100, color=[HOST_COL[g] for g in q.grp], ec=INK, lw=0.4, height=0.7, zorder=2)
    for i, r in q.iterrows():
        ax.text(r.q20_fraction * 100 + 1, i, f"{r.q20_fraction*100:.0f}%  (mean Phred {r.mean_phred:.0f})", va="center", fontsize=6.5)
    ax.set_yticks(y); ax.set_yticklabels(list(q.host), fontsize=7.5, style="italic")
    ax.invert_yaxis(); ax.set_xlim(0, 118); ax.set_xticks(range(0, 101, 20))
    ax.set_xlabel("Bases with Phred ≥ 20 in the trimmed forward read, %")
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.set_title("Sanger read quality per isolate, ranked from best to worst", loc="left", fontsize=8.5, pad=6)
    from matplotlib.patches import Patch
    used = [g for g in HOST_COL if g in set(q.grp)]
    ax.legend([Patch(color=HOST_COL[g], ec=INK, lw=0.4) for g in used], [species_label(g, groups, "sci") for g in used],
              fontsize=6.3, frameon=False, loc="upper center", bbox_to_anchor=(0.35, -0.16), ncol=3, title="Host group", title_fontsize=6.5)
    save(fig, out, "FigS2_sanger_qc")


def fig_selection(sites: pd.DataFrame, motifs: dict, out: Path, motif_table: pd.DataFrame | None = None, meta_dir: Path | None = None):
    """A: codon-wise selection tests with Rep motifs; B: residues at each motif in the ten Kerman isolates."""
    meta_dir = Path(meta_dir) if meta_dir else Path(__file__).resolve().parents[1] / "metadata"
    hosts = load_host_mapping(meta_dir)
    fig = plt.figure(figsize=(7.4, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.38)
    ax = fig.add_subplot(gs[0])
    fel = sites[sites.method == "FEL"]; me = sites[sites.method == "MEME"]
    pf = next(c for c in fel.columns if c.lower().startswith("p-value"))
    yf = -np.log10(fel[pf].astype(float).clip(lower=1e-6)); sig = fel[pf].astype(float) < 0.1
    ax.scatter(fel.codon[~sig], yf[~sig], s=9, color="#BBBBBB", label="FEL, not significant", zorder=2)
    ax.scatter(fel.codon[sig], yf[sig], s=13, color="#222222", label=f"FEL, p < 0.1 (n = {int(sig.sum())}; all purifying)", zorder=3)
    ym = -np.log10(me[pf].astype(float).clip(lower=1e-6))
    ax.scatter(me.codon, ym, s=9, marker="^", facecolor="none", edgecolor="#777777", lw=0.5, label="MEME (no site with p < 0.1)", zorder=2)
    ax.axhline(1, ls="--", color="0.5", lw=0.7); ax.text(1, 1.07, "p = 0.1", fontsize=6, color="0.4", ha="left")
    ytop = ax.get_ylim()[1]
    for i, (k, (a0, a1)) in enumerate(motifs.items()):
        ax.axvspan(a0, a1 + 0.5, color="#D9E6F2", zorder=0)
        short = k.split(" (")[0].replace("RCR motif", "RCR").replace("SF3 ", "")
        ax.text((a0 + a1) / 2, ytop * 1.02, short, fontsize=6.3, ha="center", va="bottom", fontweight="bold")
    ax.set_xlim(0, 218); ax.set_ylim(-0.1, 5.4); ytop = 5.4
    ax.set_xlabel("Codon position in the 217-codon Rep fragment"); ax.set_ylabel("−log$_{10}$ p")
    ax.legend(fontsize=6.2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, columnspacing=1.4); ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.text(-0.075, 1.12, "A", transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right")
    ax.set_title("Codon-by-codon selection (blue bands = conserved Rep motifs)", loc="left", fontsize=8.5, pad=14)
    # ---- B: motif alignment
    bx = fig.add_subplot(gs[1]); bx.axis("off")
    bx.text(-0.075, 1.07, "B", transform=bx.transAxes, fontsize=11, fontweight="bold", ha="right")
    bx.set_title("Residues at the five Rep motifs in the ten Kerman isolates (variant residues highlighted)", loc="left", fontsize=8.5, pad=6)
    if motif_table is not None and len(motif_table):
        mt = motif_table.copy()
        mcols = [c for c in mt.columns if "motif" in c.lower() or "Walker" in c]
        mt["vouch"] = mt.sequence.str.split("_Kerman_").str[-1]
        mt["sp"] = mt.vouch.map(lambda v: host_label(v, hosts))
        mt = mt.sort_values(["sp", "vouch"]).reset_index(drop=True)
        parse = lambda s: s.split("@")[0]
        cons = {c: collections.Counter(mt[c].map(parse)).most_common(1)[0][0] for c in mcols}
        n = len(mt) + 3
        colx = np.linspace(0.34, 0.97, len(mcols) + 1)[:-1]
        cw = colx[1] - colx[0]
        bx.set_xlim(0, 1); bx.set_ylim(n, -0.4)
        for j, c in enumerate(mcols):
            pos0 = int(mt[c].iloc[0].split("@")[1])
            bx.text(colx[j] + cw / 2, 0, c.replace(" (", "\n(").replace("hhhDD", "hhhDD"), ha="center", va="center", fontsize=5.8, fontweight="bold", linespacing=1.1)
            bx.text(colx[j] + cw / 2, 0.95, f"codon {pos0}", ha="center", va="center", fontsize=5.5, color="0.4")
        for i, r in enumerate(mt.itertuples()):
            y = i + 3.1
            bx.text(0.33, y, iso_label(r.vouch, hosts), ha="right", va="center", fontsize=6.3, style="italic")
            for j, c in enumerate(mcols):
                res = parse(getattr(mt.iloc[i], c)); diff = res != cons[c]
                bx.add_patch(Rectangle((colx[j] + 0.01, y - 0.4), cw - 0.02, 0.8, color="#F4B183" if diff else "#EEF3F8", ec="none"))
                bx.text(colx[j] + cw / 2, y, res, ha="center", va="center", fontsize=7, family="monospace", fontweight="bold" if diff else "normal")
        bx.text(0.33, 2.15, "consensus", ha="right", va="center", fontsize=6.3, color="0.3")
        for j, c in enumerate(mcols):
            bx.text(colx[j] + cw / 2, 2.15, cons[c], ha="center", va="center", fontsize=6.5, family="monospace", color="0.3")
    save(fig, out, "FigS3_selection_sites")
