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