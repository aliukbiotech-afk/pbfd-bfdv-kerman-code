import zipfile
if not (WORK / "inputs" / "scripts").exists():
    zipfile.ZipFile("/content/PBFD_colab_inputs.zip").extractall(WORK / "inputs")
INP = WORK / "inputs"; DATA, META = INP / "data", INP / "metadata"
sys.path.insert(0, str(INP / "scripts"))
import pbfd_phylo as P, pbfd_stats as S, pbfd_figures as F, make_figures as MF
from pbfd_names import pending_names
import numpy as np, pandas as pd
RES = WORK / "results"
for d in ["tables", "trees", "alignments", "figures", "placement", "selection", "recombination", "logs", "stats"]:
    (RES / d).mkdir(parents=True, exist_ok=True)
LOG = RES / "logs" / "commands.log"
run = lambda cmd, **kw: P.run(cmd, LOG, **kw)
print(sorted(p.name for p in DATA.iterdir())); print(pending_names(META))