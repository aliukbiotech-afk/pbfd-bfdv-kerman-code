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