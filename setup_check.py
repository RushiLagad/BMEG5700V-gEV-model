"""
One-time setup and preflight. Run from inside this folder:

    python setup_check.py

Installs what is missing, applies the PyDREAM NumPy 2 patch, then smoke-tests both
models and the saved fits. Safe to run more than once.

Nothing here needs BioNetGen or Perl. The reaction networks are pre-generated into
network_M1.py and network_M3.py, so the models run on numpy and scipy alone. PySB and
BioNetGen are only needed to regenerate those files (make_network.py) or to check them
against the PySB source (verify_network.py).
"""
import importlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

REQUIRED = [("numpy", "numpy"), ("scipy", "scipy"), ("pandas", "pandas"),
            ("matplotlib", "matplotlib"), ("pydream", "pydream")]
OPTIONAL = [("pysb", "pysb"), ("bionetgen", "bionetgen")]


def have(mod):
    try:
        importlib.import_module(mod)
        return True
    except Exception:
        return False


print("=" * 70)
print("1. packages needed to run the models")
missing = [pip for mod, pip in REQUIRED if not have(mod)]
if missing:
    print("   installing:", " ".join(missing))
    if subprocess.run([sys.executable, "-m", "pip", "install", *missing]).returncode != 0:
        print("   pip failed. If this is a Microsoft Store Python, install python.org Python "
              "or add --user, then rerun.")
        sys.exit(1)
else:
    print("   all present")
for mod, _ in REQUIRED:
    try:
        m = importlib.import_module(mod)
        print(f"   {mod:12s} {getattr(m, '__version__', 'installed')}")
    except Exception as e:
        print(f"   {mod:12s} FAILED TO IMPORT: {e}")
        sys.exit(1)

print("\n   optional, only for regenerating the networks from the PySB source:")
for mod, _ in OPTIONAL:
    print(f"   {mod:12s} {'present' if have(mod) else 'not installed, which is fine'}")

print("\n2. PyDREAM NumPy 2 patch")
import numpy as np
import pydream
dream_py = os.path.join(os.path.dirname(pydream.__file__), "Dream.py")
src = open(dream_py, encoding="utf-8").read()
bad = re.findall(r"int\(np\.where\([^\n]*?\)\[0\]\)", src)
if int(np.__version__.split(".")[0]) < 2:
    print(f"   numpy {np.__version__} is below 2, no patch needed")
elif not bad:
    print("   already patched, or this version does not need it")
else:
    open(dream_py + ".bak", "w", encoding="utf-8").write(src)
    for b in set(bad):
        src = src.replace(b, b[:-1] + "[0])")
    open(dream_py, "w", encoding="utf-8").write(src)
    print(f"   patched {len(set(bad))} line(s) in {dream_py}")
    print(f"   backup at {dream_py}.bak")

print("\n3. smoke test: both reaction networks load and integrate")
from engine import load
for name in ("M1", "M3"):
    sim = load(name)
    pv = np.tile(sim.defaults, (2, 1))
    pv[:, sim.index["gEV_0"]] = [0.0, 500.0]
    out = sim.run(pv, np.linspace(0, 600, 60))
    ratio = out[1]["obs_P"][-1] / out[0]["obs_P"][-1]
    print(f"   {name}: {sim.kind} engine, {len(sim.names)} parameters, "
          f"P at dose 500 is {ratio:.2f}x control  OK")

print("\n4. smoke test: saved fits present")
import json
for mdl in ("M1", "M3"):
    full = os.path.join(HERE, f"dream_{mdl}", "summary_full.json")
    path = full if os.path.exists(full) else os.path.join(HERE, f"dream_{mdl}", "summary.json")
    if not os.path.exists(path):
        print(f"   MISSING {path}")
        sys.exit(1)
    s = json.load(open(path))
    print(f"   dream_{mdl}: max log-lik {s['max_loglik']:.1f}  AIC {s['AIC']:.1f}  OK")

print("\n" + "=" * 70)
print("Ready. Demo order:")
print("   python demo.py")
print("   python demo.py --threshold 60")
print("   python demo.py --threshold 85")
print("   python demo_fit.py --live 150")
print("=" * 70)
