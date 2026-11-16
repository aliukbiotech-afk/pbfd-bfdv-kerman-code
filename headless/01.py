import os
#@title Config
FAST = os.environ.get("PBFD_FAST", "0") == "1"
EMAIL = ""   # blank on purpose: no personal e-mail is sent to NCBI
SEED = 20261103
import os, sys, pathlib, time, json, shutil, subprocess
IN_COLAB = True
WORK = pathlib.Path("/content/pbfd") if IN_COLAB else pathlib.Path.cwd() / "colab_work"
WORK.mkdir(parents=True, exist_ok=True); os.chdir(WORK)
THREADS = os.cpu_count() or 2
UFBOOT = 1000; ALRT = 1000; NPERM = 1000; NMC = 100_000
if FAST:
    UFBOOT, ALRT, NPERM, NMC = 1000, 0, 200, 20_000
print("work dir:", WORK, "| threads:", THREADS, "| FAST:", FAST)