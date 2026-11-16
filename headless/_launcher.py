"""Runs on the Colab VM (via `colab exec -f`). Defines run_cells(first, last) which executes the
uploaded cell scripts /content/cells/NN.py in a background thread, sharing ONE namespace (the kernel's
__main__), so a failed cell can be fixed and the run resumed with the state intact.
Progress -> /content/status.json ; per-cell output -> /content/cells_out/NN.txt
"""
import threading, traceback, json, time, io, contextlib, sys, os, pathlib

CELLS = pathlib.Path("/content/cells"); OUT = pathlib.Path("/content/cells_out"); OUT.mkdir(exist_ok=True)
STATUS = pathlib.Path("/content/status.json")


def _write(**kw):
    d = json.loads(STATUS.read_text()) if STATUS.exists() else {}
    d.update(kw, updated=time.strftime("%H:%M:%S")); STATUS.write_text(json.dumps(d))


class _Tee(io.StringIO):
    def __init__(self, path): super().__init__(); self.f = open(path, "w", buffering=1)
    def write(self, s): self.f.write(s); return super().write(s)


def _run(first, last):
    g = globals()
    for n in range(first, last + 1):
        p = CELLS / f"{n:02d}.py"
        if not p.exists():
            continue
        _write(current=n, state="running", started=time.strftime("%H:%M:%S"))
        t = _Tee(OUT / f"{n:02d}.txt")
        try:
            with contextlib.redirect_stdout(t), contextlib.redirect_stderr(t):
                exec(compile(p.read_text(), str(p), "exec"), g)
        except BaseException:
            t.f.write("\n" + traceback.format_exc()); t.f.close()
            _write(state="FAILED", failed_cell=n, error=traceback.format_exc()[-1500:])
            return
        t.f.close(); _write(last_done=n)
    _write(state="DONE")


def run_cells(first, last):
    _write(state="starting", first=first, last=last, failed_cell=None, error=None)
    th = threading.Thread(target=_run, args=(first, last), daemon=True); th.start()
    return th
