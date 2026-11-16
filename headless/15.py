mot = P.motif_table({k: aln[k].replace("-", "") for k in KER})
mot.to_csv(RES / "tables" / "rep_motifs_kerman.csv", index=False)
# subset for recombination / selection: Kerman + Iranian + 5 nearest per isolate (in-frame, no stops)
sub_names = list(dict.fromkeys(KER + IRAN + near.groupby("isolate").head(5).neighbour.tolist()))
subnt = {k: aln[k].replace("-", "") for k in sub_names}
subnt = {k: v for k, v in subnt.items() if "*" not in P.translate(v)[:-1] and len(v) >= 600}
print("selection/recombination subset:", len(subnt), "sequences")
SEL = RES / "selection"
cod = P.codon_alignment(subnt, MAFFT, SEL, LOG)
cod = {k: v[: len(v) // 3 * 3] for k, v in cod.items()}
P.write_fasta(cod, SEL / "subset_codon_aln.fasta")
P.write_fasta(cod, RES / "recombination" / "subset_codon_aln_for_RDP5.fasta")   # for optional RDP5 run on Windows
run(f"{IQTREE} -s {SEL/'subset_codon_aln.fasta'} -m GTR+F+G4 -T AUTO --threads-max {THREADS} --seed {SEED} --prefix {SEL/'subset'} -redo")
# motif codon ranges in the codon alignment (for Fig S3)
ref_aa = P.translate(cod[KER[0]].replace("-", ""))
import re as _re
cols = [i for i in range(0, len(cod[KER[0]]), 3) if cod[KER[0]][i:i+3] != "---"]
mrows = []
for k, pat in P.REP_MOTIFS.items():
    m = _re.search(pat, ref_aa)
    if m:
        mrows.append(dict(motif=k, start=cols[m.start()] // 3 + 1, end=cols[m.end() - 1] // 3 + 1))
pd.DataFrame(mrows).to_csv(RES / "tables" / "rep_motif_codons.csv", index=False)
mot