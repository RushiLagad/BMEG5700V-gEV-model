"""Simulate M1 and M3 to steady state across a dose sweep and plot the dose-response."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from engine import load

DOSES = np.array([0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500, 600])
T = np.linspace(0, 400, 400)          # long enough to reach steady state
MEASURED = [0, 150, 500]              # the three doses we actually have


def sweep(model_name):
    sim = load(model_name)
    pv = np.tile(sim.defaults, (len(DOSES), 1))
    pv[:, sim.index["gEV_0"]] = DOSES
    obs = sim.run(pv, T)
    return sim, {k: np.array([o[k][-1] for o in obs]) for k in sim.observables}


out = {}
for mname, label in (("M1", "M1_one_receptor"), ("M3", "M3_receptor_plus_brake")):
    sim, r = sweep(mname)
    out[label] = r
    print(f"\n{label}: {sim.kind} engine")
    print("  dose   occ   braked   Sig    VCAM1    P")
    for i, d in enumerate(DOSES):
        print(f"  {d:4d}  {r['obs_occ'][i]:5.1f}  {r['obs_braked'][i]:6.1f}  "
              f"{r['obs_Sig'][i]:5.1f}  {r['obs_VCAM1'][i]:6.2f}  {r['obs_P'][i]:6.2f}")

# ---------- plot ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
colors = {"M1_one_receptor": "#8A93A0", "M3_receptor_plus_brake": "#1F9E89"}
labels = {"M1_one_receptor": "M1: one receptor, no brake", "M3_receptor_plus_brake": "M3: receptor + β-arrestin2 brake"}
for name, r in out.items():
    c, ls = colors[name], "--" if name.startswith("M1") else "-"
    axes[0].plot(DOSES, r["obs_VCAM1"] / r["obs_VCAM1"][0], ls, color=c, lw=2.5, label=labels[name])
    axes[1].plot(DOSES, r["obs_P"] / r["obs_P"][0], ls, color=c, lw=2.5, label=labels[name])
    axes[2].plot(DOSES, r["obs_occ"], ls, color=c, lw=2, label=f"{name[:2]} active receptor")
    axes[2].plot(DOSES, r["obs_braked"], ":", color=c, lw=2, label=f"{name[:2]} braked receptor")
for ax, title in zip(axes[:2], ["VCAM1 (relative to control)", "Flip protein P (relative to control)"]):
    ax.axhline(1, color="#AEB7C2", lw=1, ls=":")
    for d in MEASURED:
        ax.axvline(d, color="#E8822E", lw=1, alpha=.5)
    ax.set_title(title); ax.set_xlabel("γ-EV dose (arbitrary units)"); ax.set_ylim(0, None)
axes[2].set_title("Receptor states"); axes[2].set_xlabel("γ-EV dose (arbitrary units)")
axes[0].legend(fontsize=9); axes[2].legend(fontsize=8)
fig.suptitle("Toy parameters, shape only. Orange lines = the three measured doses (0, 150, 500)", fontsize=10, color="#6B7785")
fig.tight_layout()
png = os.path.join(os.path.dirname(__file__), "dose_response_M1_vs_M3.png")
fig.savefig(png, dpi=160)
print("\nsaved", png)
