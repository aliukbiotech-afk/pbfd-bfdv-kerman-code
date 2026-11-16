#@title Cross-check: IQ-TREE with the backbone as topological constraint
run(f"{IQTREE} -s {A} -g {PL/'backbone_pruned.nwk'} -m {best_model} -B {UFBOOT} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'amplicon_constrained'} -redo")
ct = P.midpoint(P.read_tree(T / "amplicon_constrained.treefile"))
cc = P.propagate_clades(ct, lab_ref)
cmp = (aclades[aclades.tip.isin(KER)][["tip", "clade"]].rename(columns={"clade": "amplicon_tree"})
       .merge(cc[cc.tip.isin(KER)][["tip", "clade"]].rename(columns={"clade": "constrained_tree"}), on="tip")
       .merge(best[["query", "assigned_clade", "clade_LWR"]].rename(columns={"query": "tip", "assigned_clade": "EPA_ng"}), on="tip"))
cmp["agree_all"] = (cmp.amplicon_tree == cmp.constrained_tree) & (cmp.constrained_tree == cmp.EPA_ng)
cmp.to_csv(RES / "tables" / "kerman_clade_assignment_comparison.csv", index=False); cmp