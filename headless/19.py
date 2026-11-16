st = S.run_all(META, RES / "stats", NMC)
print(st["prevalence"]["overall"].to_string(index=False))
print(st["prevalence"]["by_species"][["group", "label", "n", "positive", "prevalence", "ci_low", "ci_high"]].round(3).to_string(index=False))
for k, v in st["tests"].items():
    if isinstance(v, dict):
        print(f"- {k}:", {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "table"})
st["regression"].round(3)