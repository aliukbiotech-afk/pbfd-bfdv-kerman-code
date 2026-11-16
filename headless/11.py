#@title EPA-ng placement of the 10 partial sequences on the full-genome backbone
PL = RES / "placement"
refs_bb = [t.name for t in bb.get_terminals() if t.name in aln]
P.write_fasta({k: aln[k] for k in refs_bb}, PL / "ref_amplicon_msa.fasta")
P.write_fasta({k: aln[k] for k in KER}, PL / "query.fasta")
import re
nwk = open(T / "backbone_fullgenome_MFP.treefile").read()
nwk = re.sub(r"\)[0-9.\/]+:", "):", nwk)                       # strip support labels
from Bio import Phylo
tt = Phylo.read(__import__("io").StringIO(nwk), "newick")
for tip in [x for x in tt.get_terminals() if x.name not in set(refs_bb)]:
    tt.prune(tip)
Phylo.write(tt, PL / "backbone_pruned.nwk", "newick")
(PL / "backbone_pruned.nwk").write_text(re.sub(r"\)[0-9.]+:", "):", (PL / "backbone_pruned.nwk").read_text()))
run(f"{RAXML} --evaluate --msa {PL/'ref_amplicon_msa.fasta'} --tree {PL/'backbone_pruned.nwk'} --model GTR+G --prefix {PL/'eval'} --threads {THREADS} --redo --seed {SEED}")
run(f"{EPA} --ref-msa {PL/'ref_amplicon_msa.fasta'} --tree {PL/'eval.raxml.bestTree'} --query {PL/'query.fasta'} "
    f"--model {PL/'eval.raxml.bestModel'} --outdir {PL} --redo --threads {THREADS}")
ptree, pl = P.parse_jplace(PL / "epa_result.jplace")
tipc = dict(zip(bclades.tip, bclades.clade))
pl, summ, best = P.placement_clade_summary(ptree, pl, tipc)
pl.to_csv(RES / "tables" / "placement_all.csv", index=False)
summ.to_csv(RES / "tables" / "placement_lwr_by_clade.csv", index=False)
best.to_csv(RES / "tables" / "placement_best.csv", index=False)
best