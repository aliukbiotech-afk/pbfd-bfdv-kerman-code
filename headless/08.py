A = RES / "alignments" / "amplicon_all_aln.fasta"; T = RES / "trees"
alrt = f"-alrt {ALRT}" if ALRT else ""
run(f"{IQTREE} -s {A} -m MFP -B {UFBOOT} {alrt} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'amplicon_MFP'} -redo")
run(f"{IQTREE} -s {A} -m TIM3+F+I+G4 -B {UFBOOT} -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {T/'amplicon_TIM3'} -redo")
best_model = [l for l in open(T / "amplicon_MFP.iqtree") if l.startswith("Best-fit model")][0].split(":")[1].strip()
run(f"{IQTREE} -rf {T/'amplicon_MFP.treefile'} {T/'amplicon_TIM3.treefile'} --prefix {T/'rf_MFP_vs_TIM3'} -redo")
print("best-fit model (BIC):", best_model)
print(open(T / "rf_MFP_vs_TIM3.rfdist").read())