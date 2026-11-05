"""
Diagnose and fix the BioNetGen / perl problem on Windows. Run:

    python fix_bng.py

It reports which perl is being used, finds the bundled BioNetGen, then tries to
generate the reaction network under several configurations and tells you which
one worked. Nothing is installed or changed on disk.
"""
import os
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

print("=" * 70)
print("1. perl")
perl = shutil.which("perl")
print("   which perl:", perl or "NOT FOUND ON PATH")
if perl:
    try:
        out = subprocess.run([perl, "-e", "print $^O, '|', $], '|', join(';', @INC)"],
                             capture_output=True, text=True, timeout=30).stdout
        osname, ver, inc = out.split("|", 2)
        print(f"   $^O = {osname}   version {ver}")
        print("   @INC first entries:", "; ".join(inc.split(";")[:3]))
        if sys.platform == "win32" and osname != "MSWin32":
            print("   ^ unix-style perl (MSYS / Git / Cygwin) on Windows. That is the usual cause.")
    except Exception as e:
        print("   could not query perl:", e)
for cand in (r"C:\Strawberry\perl\bin\perl.exe", r"C:\Perl64\bin\perl.exe"):
    if pathlib.Path(cand).exists():
        print("   also found a native perl at:", cand)

print("\n2. bundled BioNetGen")
try:
    import bionetgen
    root = pathlib.Path(bionetgen.__file__).parent
    print("   package:", root)
    for name in ("bng-win", "bng-linux", "bng-mac"):
        d = root / name
        if d.exists():
            ok = (d / "Perl2" / "BNGUtils.pm").exists()
            print(f"   {name:10s} exists, Perl2/BNGUtils.pm {'present' if ok else 'MISSING'}")
except ImportError:
    print("   bionetgen not installed"); sys.exit(1)

print("\n3. trying to generate the network")


def attempt(label, env_updates):
    saved = {k: os.environ.get(k) for k in env_updates}
    os.environ.update({k: v for k, v in env_updates.items() if v is not None})
    for k, v in env_updates.items():
        if v is None:
            os.environ.pop(k, None)
    # a fresh interpreter each time, so a half-initialised pysb never carries over
    code = (
        "import sys; sys.path.insert(0, r'%s');"
        "from gev_models import build;"
        "from pysb.bng import generate_equations;"
        "m = build(with_brake=True); generate_equations(m);"
        "print('SPECIES', len(m.species), 'REACTIONS', len(m.reactions))" % str(HERE)
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600)
    for k, v in saved.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    ok = "SPECIES" in r.stdout
    print(f"   {label:38s} {'OK  ' + r.stdout.strip().splitlines()[-1] if ok else 'failed'}")
    if not ok:
        tail = [l for l in (r.stderr or "").strip().splitlines() if l.strip()][-2:]
        for l in tail:
            print("       ", l[:160])
    return ok


winpath = root / ("bng-win" if sys.platform == "win32" else ("bng-mac" if sys.platform == "darwin" else "bng-linux"))
results = {}
results["as-is"] = attempt("no BNGPATH (current behaviour)", {"BNGPATH": None})
if not results["as-is"]:
    results["fwd"] = attempt("BNGPATH with forward slashes", {"BNGPATH": winpath.as_posix()})
    if not results.get("fwd"):
        results["back"] = attempt("BNGPATH with backslashes", {"BNGPATH": str(winpath)})

print("\n" + "=" * 70)
if results.get("as-is"):
    print("BioNetGen already works. The earlier failure was something else; paste this output.")
elif results.get("fwd"):
    print("FIXED by BNGPATH with forward slashes. bng_setup.py already does this,")
    print("and demo.py / calibrate.py / run_dose_response.py import it, so just rerun:")
    print("    python demo.py")
elif results.get("back"):
    print("FIXED by BNGPATH with backslashes. Tell Claude: the posix form did not work.")
else:
    print("Still failing. Paste this whole output. The likely next step is a native perl:")
    print("    conda install -c conda-forge strawberryperl")
print("=" * 70)
