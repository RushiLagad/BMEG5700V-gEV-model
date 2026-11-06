"""
Live demo of the calibrated models.

    python demo_fit.py                 # load the saved posteriors, print the comparison, draw the figure (seconds)
    python demo_fit.py --live 150      # first run a short PyDREAM burst live (about a minute), then the above
    python demo_fit.py --no-plot

Needs dream_M1/ and dream_M3/ from a full calibrate.py run in the same folder.
"""
import argparse, json, os, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

ap = argparse.ArgumentParser()
ap.add_argument("--live", type=int, default=0, help="run this many PyDREAM iterations live for M3 before loading the saved fit")
ap.add_argument("--no-plot", action="store_true")
a = ap.parse_args()

if a.live:
    print(f"\n=== live: PyDREAM on M3, 3 chains x {a.live} iterations (watch the log-likelihood climb) ===")
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, "calibrate.py"), "--model", "M3", "--niter", str(a.live), "--nchains", "3", "--quick", "--seed", "7"],
                       capture_output=True, text=True)
    lines = [l for l in r.stdout.splitlines() if l.strip() and not l.startswith("Saving")]
    print("\n".join(lines[:3])); print("   ..."); print("\n".join(lines[-6:]))
    print(f"(that took {time.time()-t0:.0f} s; the saved fit below used 1500 iterations)")

import calibrate as C


def ref(model, stem, ext):
    """Prefer <stem>_full.<ext>, the reference run, over whatever a later run left behind."""
    full = os.path.join(HERE, f"dream_{model}", f"{stem}_full.{ext}")
    return full if os.path.exists(full) else os.path.join(HERE, f"dream_{model}", f"{stem}.{ext}")


S = {m: json.load(open(ref(m, "summary", "json"))) for m in ("M1", "M3")}

print("\n=== saved fit: M1 vs M3 against the real log2 fold changes ===")
keys = ["VCAM1@150", "VCAM1@500", "P@150", "P@500"]
print(f"{'':22s}" + "".join(f"{k:>12s}" for k in keys) + f"{'max log-lik':>14s}{'AIC':>8s}")
print(f"{'data':22s}" + "".join(f"{S['M3']['data'][k][0]:12.2f}" for k in keys))
for m, label in (("M1", "M1 one receptor"), ("M3", "M3 receptor + brake")):
    print(f"{label:22s}" + "".join(f"{S[m]['best_pred'][k]:12.2f}" for k in keys) + f"{S[m]['max_loglik']:14.1f}{S[m]['AIC']:8.1f}")
print(f"\nDelta AIC = {S['M1']['AIC'] - S['M3']['AIC']:.0f} in favour of M3, after paying for {S['M3']['n_params'] - S['M1']['n_params']} extra parameters.")
print("M1's best effort is a flat P: it cannot go down at 150 and up at 500 at ANY parameter values.")

print("\n=== what the three doses pin down (M3, 95% posterior, log10) ===")
for n, q in S["M3"]["posterior_quantiles_log10"].items():
    w = q[2] - q[0]; tag = "PINNED" if w < 1.0 else ("loose" if w < 2.5 else "unidentified")
    print(f"  {n:9s} [{q[0]:6.2f}, {q[2]:6.2f}]   {tag}")
gr = S["M3"].get("gelman_rubin")
if isinstance(gr, dict):
    print(f"\nGelman-Rubin (want < 1.1): max = {max(gr.values()):.2f}  ->  " + ("converged" if max(gr.values()) < 1.1 else "NOT yet converged; this is a 1500-iteration run, 5000+ needed before quoting intervals"))

if a.no_plot: sys.exit(0)
import matplotlib
if not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
elif matplotlib.get_backend().lower() in ("agg", "template"):
    for _bk in ("TkAgg", "QtAgg", "Qt5Agg"):
        try:
            matplotlib.use(_bk); break
        except Exception:
            continue
import matplotlib.pyplot as plt
DOSES = np.array([0,25,50,75,100,150,200,250,300,350,400,450,500,550,600], float)
rng = np.random.default_rng(0)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for model, color, label in (("M1", "#8A93A0", "M1: one receptor"), ("M3", "#1F9E89", "M3: receptor + β-arrestin2 brake")):
    C.setup(model)
    smp = np.load(ref(model, "samples_log10", "npy")); lp = np.load(ref(model, "logp", "npy"))
    def curves(theta):
        pv = np.tile(C._S["base"], (len(DOSES), 1)); pv[:, C._S["idx"]] = 10 ** theta; pv[:, C._S["gidx"]] = DOSES
        obs = C._S["sim"].run(pv, C.T)
        V = np.array([o["obs_VCAM1"][-1] for o in obs]); P = np.array([o["obs_P"][-1] for o in obs])
        return np.log2(V / V[0]), np.log2(P / P[0])
    for th in smp[rng.choice(len(smp), 40, replace=False)]:
        v, p = curves(th); ax[0].plot(DOSES, v, color=color, alpha=0.12, lw=1.2); ax[1].plot(DOSES, p, color=color, alpha=0.12, lw=1.2)
    v, p = curves(smp[np.argmax(lp)]); ax[0].plot(DOSES, v, color=color, lw=3, label=label); ax[1].plot(DOSES, p, color=color, lw=3, label=label)
for axx, key, title in ((ax[0], "VCAM1", "VCAM1, log2 fold change vs control"), (ax[1], "P", "flip proteins (mean of 5 anchors), log2FC vs control")):
    for d in (150, 500):
        mu, sd = C.DATA[(key, d)]; axx.errorbar([d], [mu], yerr=[sd], fmt="o", color="#16222E", ms=8, capsize=5, zorder=5, label="data" if d == 150 else None)
    axx.axhline(0, color="#AEB7C2", ls=":"); axx.set_title(title, fontsize=11); axx.set_xlabel("γ-EV dose (mg/kg/day)")
ax[1].legend(fontsize=9, loc="upper left")
fig.suptitle("Thin lines: posterior draws. Thick: best fit. Black: data. Error bars are a placeholder sigma until limma SEs are exported.", fontsize=9, color="#6B7785")
fig.tight_layout(); out = os.path.join(HERE, "demo_fit.png"); fig.savefig(out, dpi=150); print("\nfigure saved:", out, f"  [matplotlib backend: {matplotlib.get_backend()}]")
try: plt.show()
except Exception: pass
