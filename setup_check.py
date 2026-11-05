"""
One-time setup and preflight for the gamma-EV model demo. Run from inside this folder:

    python setup_check.py

Installs anything missing, applies the PyDREAM NumPy 2 patch, then runs a smoke test
of both models and the saved fit. Safe to run more than once.
"""
import importlib, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKGS = [("pysb", "pysb"), ("bionetgen", "bionetgen"), ("numpy", "numpy"),
        ("pandas", "pandas"), ("matplotlib", "matplotlib"), ("scipy", "scipy"), ("pydream", "pydream")]

def have(mod):
    try:
        importlib.import_module(mod); return True
    except Exception:
        return False

print("=" * 70)
print("1. packages")
missing = [pip for mod, pip in PKGS if not have(mod)]
if missing:
    print("   installing:", " ".join(missing))
    r = subprocess.run([sys.executable, "-m", "pip", "install", *missing])
    if r.returncode != 0:
        print("   pip failed. Try: conda activate <your env>  then rerun."); sys.exit(1)
else:
    print("   all present")
for mod, _ in PKGS:
    try:
        m = importlib.import_module(mod)
        print(f"   {mod:12s} {getattr(m, '__version__', '?')}")
    except Exception as e:
        print(f"   {mod:12s} FAILED TO IMPORT: {e}"); sys.exit(1)

print("\n2. PyDREAM NumPy 2 patch")
import numpy as np
import pydream
dream_py = os.path.join(os.path.dirname(pydream.__file__), "Dream.py")
src = open(dream_py, encoding="utf-8").read()
# the two offending lines call int() on a 1-element array, which NumPy 2 refuses
bad = re.findall(r"int\(np\.where\([^\n]*?\)\[0\]\)", src)
if int(np.__version__.split(".")[0]) < 2:
    print(f"   numpy {np.__version__} < 2, no patch needed")
elif not bad:
    print("   already patched (or this version does not need it)")
else:
    open(dream_py + ".bak", "w", encoding="utf-8").write(src)
    for b in set(bad):
        src = src.replace(b, b[:-1] + "[0])")
    open(dream_py, "w", encoding="utf-8").write(src)
    print(f"   patched {len(set(bad))} line(s) in {dream_py}")
    print(f"   backup at {dream_py}.bak")

print("\n3. smoke test: build both models")
sys.path.insert(0, HERE)
from gev_models import build
for name, brake in (("M1", False), ("M3", True)):
    m = build(with_brake=brake)
    print(f"   {name}: {len(m.rules)} rules, {len(m.parameters)} parameters  OK")

print("\n4. smoke test: saved fit present")
import json
for mdl in ("M1", "M3"):
    p = os.path.join(HERE, f"dream_{mdl}", "summary.json")
    if not os.path.exists(p):
        print(f"   MISSING {p}"); sys.exit(1)
    s = json.load(open(p))
    print(f"   dream_{mdl}/summary.json  max log-lik {s['max_loglik']:.1f}  AIC {s['AIC']:.1f}  OK")

print("\n" + "=" * 70)
print("Ready. Demo order:")
print("   python demo.py")
print("   python demo.py --threshold 60")
print("   python demo.py --threshold 85")
print("   python demo_fit.py --live 150")
print("=" * 70)
