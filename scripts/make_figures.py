"""Re-render every figure from saved results + metadata (no heavy computation).

Use after filling FINAL_HOST / FINAL_SPECIES in analysis/metadata/*.csv:
    python analysis/scripts/make_figures.py  [results_dir]  [metadata_dir]
The Colab notebook calls render_all() too, so figures are identical in both places.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pbfd_stats as S
import pbfd_figures as F
import pbfd_phylo as P
from pbfd_names import pending_names, species_legend_table, load_host_mapping


def _exists(*paths):
    return all(Path(p).exists() for p in paths)


def render_all(res: Path, meta: Path) -> list[str]:
    res, meta = Path(res), Path(meta)
    fig = res / "figures"; tab = res / "tables"; done = []
    # Fig 1-2: epidemiology (fast, recomputed from metadata so new names propagate)
    birds = S.load_birds(meta); prev = S.prevalence_tables(birds)
    F.fig_prevalence(prev, meta, fig); F.fig_alluvial(birds, meta, fig); done += ["Fig1", "Fig2"]
    refmd = {}
    if (tab / "reference_metadata_merged.csv").exists():
        md = pd.read_csv(tab / "reference_metadata_merged.csv")
        refmd = {r["accession"]: r for r in md.to_dict("records")}
    # Fig 3 / S1: ML trees
    tc = dict(pd.read_csv(tab / "tip_clades_amplicon.csv")[["tip", "clade"]].values) if (tab / "tip_clades_amplicon.csv").exists() else {}
    tf = res / "trees" / "amplicon_MFP.treefile"
    if tf.exists():
        tree = P.midpoint(P.read_tree(tf))
        keep = [l.strip() for l in open(tab / "focus_tips.txt")] if (tab / "focus_tips.txt").exists() else None
        if keep:
            F.fig_tree_focus(tree, keep, tc, refmd, meta, fig); done.append("Fig3")
        F.fig_tree_full(P.midpoint(P.read_tree(tf)), tc, refmd, meta, fig); done.append("FigS1")
    bt = res / "trees" / "backbone_fullgenome_MFP.treefile"
    if bt.exists() and (tab / "tip_clades_backbone.csv").exists():
        bc = dict(pd.read_csv(tab / "tip_clades_backbone.csv")[["tip", "clade"]].values)
        F.fig_tree_full(P.midpoint(P.read_tree(bt)), bc, refmd, meta, fig, name="FigS4_backbone_fullgenome_tree")
        done.append("FigS4")
    # Fig 4: placement
    if (tab / "placement_lwr_by_clade.csv").exists():
        F.fig_placement(pd.read_csv(tab / "placement_lwr_by_clade.csv"), meta, fig); done.append("Fig4")
    # Fig 5: identity + network
    if _exists(tab / "identity_matrix_kerman_iran.csv", res / "alignments" / "amplicon_all_aln.fasta", tab / "network_members.csv"):
        D = 1 - pd.read_csv(tab / "identity_matrix_kerman_iran.csv", index_col=0) / 100
        aln = P.read_fasta(res / "alignments" / "amplicon_all_aln.fasta")
        mem = pd.read_csv(tab / "network_members.csv")
        hap, _ = P.haplotypes(aln, list(mem.tip))
        G = P.haplotype_msn(hap["seq_by_hap"])
        F.fig_identity_network(D, G, hap["hap_of"], hap["counts"], dict(mem[["tip", "group"]].values), meta, fig)
        done.append("Fig5")
    # supplementary
    if (tab / "sanger_qc.csv").exists():
        F.fig_sanger_qc(pd.read_csv(tab / "sanger_qc.csv"), fig); done.append("FigS2")
    if (tab / "selection_sites.csv").exists():
        mot = pd.read_csv(tab / "rep_motif_codons.csv") if (tab / "rep_motif_codons.csv").exists() else pd.DataFrame(columns=["motif", "start", "end"])
        F.fig_selection(pd.read_csv(tab / "selection_sites.csv"), {r.motif: (r.start, r.end) for r in mot.itertuples()}, fig)
        done.append("FigS3")
    species_legend_table(meta).to_csv(fig / "LEGEND_placeholder_names.csv", index=False)
    load_host_mapping(meta).to_csv(fig / "LEGEND_isolate_hosts.csv", index=False)
    print("rendered:", done, "| pending names:", pending_names(meta))
    return done


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    res = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "results"
    meta = Path(sys.argv[2]) if len(sys.argv) > 2 else root / "metadata"
    render_all(res, meta)
