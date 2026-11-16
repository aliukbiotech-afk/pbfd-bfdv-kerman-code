#@title Software versions
vers = {}
for name, cmd in [("mafft", f"{MAFFT} --version"), ("iqtree", f"{IQTREE} --version"), ("raxml-ng", f"{RAXML} --version"),
                  ("epa-ng", f"{EPA} --version"), ("hyphy", f"{HYPHY} --version")]:
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    vers[name] = (out.stdout + out.stderr).strip().splitlines()[0:3]
import Bio, scipy, statsmodels, matplotlib, networkx
vers.update(python=sys.version.split()[0], biopython=Bio.__version__, numpy=np.__version__, pandas=pd.__version__,
            scipy=scipy.__version__, statsmodels=statsmodels.__version__, matplotlib=matplotlib.__version__, networkx=networkx.__version__,
            seed=SEED, FAST=FAST, ufboot=UFBOOT, alrt=ALRT, permutations=NPERM, monte_carlo=NMC)
json.dump(vers, open(RES / "versions.json", "w"), indent=2)
(RES / "versions.txt").write_text("\n".join(f"{k}: {v}" for k, v in vers.items()))
print(open(RES / "versions.txt").read())