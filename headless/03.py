#@title Get inputs (cloned from GitHub)
REPO_URL = "https://github.com/aliukbiotech-afk/pbfd-bfdv-kerman-code"
INP = WORK / "inputs"
if not IN_COLAB and (WORK.parent / "scripts" / "pbfd_phylo.py").exists():
    INP = WORK.parent                      # notebook run inside a local clone of the repository
elif not (INP / "scripts").exists():
    subprocess.run(["git", "clone", "-q", "--depth", "1", REPO_URL, str(INP)], check=True)
DATA, META = INP / "data", INP / "metadata"
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