summary = dict(n_birds=136, n_positive=int(st["birds"].positive.sum()), n_isolates=len(KER),
               amplicon_alignment=[len(aln), len(next(iter(aln.values())))], best_model=best_model,
               placement=best[["query", "assigned_clade", "clade_LWR"]].to_dict("records"),
               pending_names=pending_names(META))
json.dump(summary, open(RES / "RUN_SUMMARY.json", "w"), indent=2, default=str)
shutil.make_archive(str(WORK / "results"), "zip", RES)
print(json.dumps(summary, indent=1, default=str)[:3000])
shutil.copy(WORK / "results.zip", "/content/results.zip"); print("results.zip ->", "/content/results.zip")