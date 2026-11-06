"""
Live demo: M1 (one receptor) vs M3 (receptor + beta-arrestin2 brake).

    python demo.py                      # default brake threshold (Km_arr = 72 occupied receptors)
    python demo.py --threshold 60       # move the brake; watch the flip dose move
    python demo.py --threshold 85 --hill 3
    python demo.py --no-plot            # table only

Runs both models across 13 doses in one batched ODE call each (a few seconds).
"""
import argparse, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from engine import load

ap = argparse.ArgumentParser()
ap.add_argument("--threshold", type=float, default=72.0, help="occupied-receptor level at half-max braking (Km_arr)")
ap.add_argument("--hill", type=float, default=12.0, help="steepness of the brake switch (n_hill)")
ap.add_argument("--no-plot", action="store_true")
args = ap.parse_args()

DOSES = np.array([0, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 500, 600], dtype=float)
T = np.linspace(0, 400, 80)

t0 = time.time()
results = {}
engine_kind = None
for label, mname in (("M1_one_receptor", "M1"), ("M3_receptor_plus_brake", "M3")):
    sim = load(mname)
    engine_kind = sim.kind
    pv = np.tile(sim.defaults, (len(DOSES), 1))
    pv[:, sim.index["gEV_0"]] = DOSES
    pv[:, sim.index["Km_arr"]] = args.threshold
    pv[:, sim.index["n_hill"]] = args.hill
    obs = sim.run(pv, T)
    results[label] = {k: np.array([o[k][-1] for o in obs]) for k in sim.observables}
print(f"built + simulated 2 models x {len(DOSES)} doses in {time.time()-t0:.1f} s [{engine_kind} engine]\n")

# ---------- table at the three measured doses ----------
def rel(r, key):
    return r[key] / r[key][0]

idx = [list(DOSES).index(d) for d in (0.0, 150.0, 500.0)]
print(f"{'':28s}{'dose 0':>10s}{'dose 150':>10s}{'dose 500':>10s}   reads as")
for name, label in (("M1_one_receptor", "M1 one receptor"), ("M3_receptor_plus_brake", "M3 receptor + brake")):
    r = results[name]
    v = rel(r, "obs_VCAM1"); p = rel(r, "obs_P")
    print(f"{label:28s}")
    print(f"{'   VCAM1 (x control)':28s}" + "".join(f"{v[i]:10.2f}" for i in idx) + "   " + ("down, down" if v[idx[2]] < 1 else "up at high"))
    shape = "down, then UP" if (p[idx[1]] < 1 and p[idx[2]] > 1) else ("down, down" if p[idx[2]] < 1 else "other")
    print(f"{'   flip protein P (x control)':28s}" + "".join(f"{p[i]:10.2f}" for i in idx) + f"   {shape}")
print(f"\nData: VCAM1 ~0.4x at 150 and down at 500 (western blot); the 8 flip proteins ~0.3x at 150 and 2 to 4x at 500.")
p3 = rel(results["M3_receptor_plus_brake"], "obs_P")
cross = None
for i in range(1, len(DOSES)):                      # linear interpolation, not grid snapping
    if p3[i - 1] < 1 <= p3[i]:
        f = (1 - p3[i - 1]) / (p3[i] - p3[i - 1])
        cross = DOSES[i - 1] + f * (DOSES[i] - DOSES[i - 1])
        break
print(f"M3 predicts P crosses back above control near dose {cross:.0f} (brake Km = {args.threshold:.0f} occupied receptors)." if cross else "M3 does not cross control in this range.")

if args.no_plot:
    sys.exit(0)

import matplotlib
if not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
elif matplotlib.get_backend().lower() in ("agg", "template"):
    for _bk in ("TkAgg", "QtAgg", "Qt5Agg"):   # pick a window backend if one is installed
        try:
            matplotlib.use(_bk); break
        except Exception:
            continue
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for name, c, ls, lab in (("M1_one_receptor", "#8A93A0", "--", "M1: one receptor"),
                         ("M3_receptor_plus_brake", "#1F9E89", "-", "M3: receptor + brake")):
    r = results[name]
    ax[0].plot(DOSES, rel(r, "obs_VCAM1"), ls, color=c, lw=2.6, label=lab)
    ax[1].plot(DOSES, rel(r, "obs_P"), ls, color=c, lw=2.6, label=lab)
for a, t in zip(ax, ("VCAM1, relative to control", "flip protein P, relative to control")):
    a.axhline(1, color="#AEB7C2", ls=":", lw=1)
    for d in (150, 500):
        a.axvline(d, color="#E8822E", lw=1, alpha=.6)
    a.set_title(t); a.set_xlabel("γ-EV dose (arbitrary units)"); a.set_ylim(0, None)
ax[0].legend()
fig.suptitle(f"Toy parameters, shape only. Orange = measured doses 150 and 500. Brake Km = {args.threshold:.0f} occupied receptors (of 100).", fontsize=10, color="#6B7785")
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"demo_threshold{int(args.threshold)}.png")
fig.savefig(out, dpi=150)
print("figure saved:", out, f"  [matplotlib backend: {matplotlib.get_backend()}]")
try:
    plt.show()
except Exception:
    pass
