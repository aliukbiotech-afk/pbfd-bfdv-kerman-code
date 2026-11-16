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