"""Split PBFD_analysis.ipynb into one script per code cell for `colab exec` (headless runs).
Interactive steps (file download, display, e-mail prompt) are replaced.
    python make_headless.py
"""
import json, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
nb = json.load(open(HERE / "PBFD_analysis.ipynb", encoding="utf-8"))
out = HERE / "headless"; out.mkdir(exist_ok=True)
for f in out.glob("[0-9][0-9].py"): f.unlink()
n = 0
for c in nb["cells"]:
    if c["cell_type"] != "code": continue
    src = "".join(c["source"]); n += 1
    first = src.splitlines()[0]
    if first.startswith("#@title Config"):
        src = src.replace('IN_COLAB = "google.colab" in sys.modules', "IN_COLAB = True")
        src = src.replace('FAST = False          #@param {type:"boolean"}', 'FAST = os.environ.get("PBFD_FAST", "0") == "1"')
        src = src.replace("import os, sys,", "import os, sys,") 
        src = src.replace('EMAIL = ""            #@param {type:"string"}', 'EMAIL = ""   # blank on purpose: no personal e-mail is sent to NCBI')
        src = "import os\n" + src
    elif 'assert EMAIL' in src:
        src = re.sub(r'assert EMAIL.*\n', '', src)
    elif "MF.render_all(RES, META)" in src:
        src = src.replace("from IPython.display import Image, display\n", "").replace(
            'for f in sorted((RES / "figures").glob("Fig*.png")):\n    print(f.name); display(Image(str(f), width=900))',
            'print(sorted(p.name for p in (RES / "figures").glob("Fig*.png")))')
    elif "files.download" in src:
        src = src.replace('if IN_COLAB:\n    from google.colab import files; files.download(str(WORK / "results.zip"))',
                          'shutil.copy(WORK / "results.zip", "/content/results.zip"); print("results.zip ->", "/content/results.zip")')
    (out / f"{n:02d}.py").write_text(src, encoding="utf-8")
print(n, "cell scripts written to", out)
