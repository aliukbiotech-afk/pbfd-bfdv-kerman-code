#@title Clade labels on the amplicon tree (same propagation, for comparison)
at = P.midpoint(P.read_tree(T / "amplicon_MFP.treefile"))
lab_all = dict(zip(bclades.tip, bclades.clade))
lab_ref = {k: v for k, v in lab_all.items() if v != "unassigned"}
aclades = P.propagate_clades(at, lab_ref)
aclades.to_csv(RES / "tables" / "tip_clades_amplicon.csv", index=False)
aclades[aclades.tip.isin(KER)]