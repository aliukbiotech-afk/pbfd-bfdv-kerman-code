# Run on the Colab VM with `colab exec -f headless/_start_blocking.py`: clones this repository, copies the
# cell scripts headless/NN.py to /content/cells/ and runs cells 1-21 in the main kernel thread.
import subprocess, shutil, pathlib
REPO = pathlib.Path("/content/pbfd-bfdv-kerman-code")
if not REPO.exists():
    subprocess.run(["git", "clone", "-q", "--depth", "1", "https://github.com/aliukbiotech-afk/pbfd-bfdv-kerman-code", str(REPO)], check=True)
pathlib.Path("/content/cells").mkdir(exist_ok=True)
for p in sorted((REPO / "headless").glob("[0-9][0-9].py")):
    shutil.copy(p, pathlib.Path("/content/cells") / p.name)
exec(open(REPO / "headless" / "_launcher.py").read())
# run in the MAIN kernel thread so Colab sees a continuously executing command
_write(state="starting", first=1, last=21, failed_cell=None, error=None)
_run(1, 21)
