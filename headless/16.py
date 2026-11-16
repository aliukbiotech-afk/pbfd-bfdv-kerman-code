#@title GARD recombination screen (HyPhy)
REC = RES / "recombination"
try:
    run(f"{HYPHY} gard --alignment {SEL/'subset_codon_aln.fasta'} --type nucleotide --rv GDD --rate-classes 3 "
        f"--output {REC/'gard.json'} --output-lf {REC/'gard_best.nex'} CPU={THREADS}")
    G = json.load(open(REC / "gard.json"))
    bp = G.get("breakpointData", {})
    print("GARD: potential breakpoints:", G.get("potentialBreakpoints"), "| best model breakpoints:", len(bp) - 1 if bp else 0)
    print(json.dumps({k: G[k] for k in G if k in ("improvements", "siteBreakPointSupport")}, indent=1)[:3000])
except Exception as e:
    print("GARD failed:", e)