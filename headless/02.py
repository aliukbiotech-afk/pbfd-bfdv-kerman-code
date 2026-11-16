#@title Install tools (micromamba + bioconda; ~5 min)
ENV = pathlib.Path("/content/env") if IN_COLAB else pathlib.Path.home() / "pbfd_env"
BIN = ENV / "bin"
if not (BIN / "iqtree2").exists() and not (BIN / "iqtree").exists():
    subprocess.run("curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba", shell=True, check=True)
    subprocess.run(f"micromamba create -y -q -p {ENV} -c conda-forge -c bioconda "
                   "mafft iqtree raxml-ng epa-ng hyphy", shell=True, check=True)
subprocess.run("pip -q install biopython networkx", shell=True, check=True)
IQTREE = str(BIN / "iqtree2") if (BIN / "iqtree2").exists() else str(BIN / "iqtree")
MAFFT, RAXML, EPA, HYPHY = (str(BIN / x) for x in ("mafft", "raxml-ng", "epa-ng", "hyphy"))