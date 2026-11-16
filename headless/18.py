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