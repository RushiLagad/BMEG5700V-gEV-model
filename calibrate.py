"""
PyDREAM calibration of M1 and M3 against the real log2 fold changes in targets.csv.

    python calibrate.py --model M3 --niter 3000 --nchains 3
    python calibrate.py --model M1 --niter 3000 --nchains 3
    python calibrate.py --model M3 --niter 100 --nchains 3 --quick     # smoke test

What is fitted
  * Steady-state readouts at dose 0, 150, 500 (one time point in the data, so steady state).
  * VCAM1 log2FC at 150 and 500.
  * P log2FC at 150 and 500, where the data value is the MEAN log2FC of the anchor flip
    proteins marked use_in_fit=yes (Nceh1, Comp, Serinc3, Aldh1l2, Sh3bgr) and the error
    is their spread. Tnnc1 is excluded pending histology; Saa1/2 are not artery-wall readouts.
  * Likelihood: Gaussian on the log2 scale. sigma = max(sd across anchors, SIGMA_FLOOR);
    SIGMA_FLOOR is a placeholder until limma's per-protein standard errors are exported.

Parameters fitted (log10 scale, uniform priors over the ranges below); all others fixed.
Dose enters as gEV_0 in the same arbitrary units as kr/kf, so the receptor Kd is fitted too.
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import load

HERE = os.path.dirname(os.path.abspath(__file__))
DOSES = np.array([0.0, 150.0, 500.0])
T = np.linspace(0, 600, 60)
SIGMA_FLOOR = 0.35      # log2 units; replace with limma SEs when available

# ---------------- data ----------------
tg = pd.read_csv(os.path.join(HERE, "targets.csv"))
tg = tg[tg.use_in_fit == "yes"]
def target(role, dose):
    rows = tg[(tg.role == role) & (tg.dose_mg_kg == dose)].log2FC_vs_control.values
    return float(rows.mean()), float(max(rows.std(ddof=0) if len(rows) > 1 else 0.0, SIGMA_FLOOR))
DATA = {  # (mean log2FC, sigma)
    ("VCAM1", 150): target("readout_vcam1", 150), ("VCAM1", 500): target("readout_vcam1", 500),
    ("P", 150): target("readout_P", 150),         ("P", 500): target("readout_P", 500),
}

# ---------------- which parameters to fit, and their log10 prior ranges ----------------
FIT = {                      # name: (lo, hi) in log10
    "kr":       (0.5, 3.5),   # Kd = kr/kf with kf fixed at 1 -> Kd in [3, 3000] dose units
    "ksig":     (-1, 1.5),
    "KI_v":     (-0.5, 2.5),
    "KI_p":     (-0.5, 2.5),
    "fold_p":   (0, 1.5),
    "KA_p":     (-0.5, 2.5),
}
FIT_M3_EXTRA = {
    "kcat_arr": (-2, 1.5),
    "Km_arr":   (0.5, 2.0),   # 3 to 100 occupied receptors (of 100)
    "n_hill":   (0, 1.2),     # 1 to ~16
    "kreset":   (-1.5, 1),
}

# ---------------- module-level state (PyDREAM pickles the likelihood, so no closures) ----------------
_S = {}

def setup(model_name):
    sim = load(model_name)
    fit = dict(FIT); fit.update(FIT_M3_EXTRA if model_name == "M3" else {})
    fit_names = list(fit)
    _S.update(sim=sim, base=sim.defaults.copy(), fit=fit, fit_names=fit_names,
              idx=[sim.index[n] for n in fit_names], gidx=sim.index["gEV_0"])
    return sim, fit_names, fit


def simulate(theta_log10):
    pv = np.tile(_S["base"], (3, 1))
    pv[:, _S["idx"]] = 10 ** np.asarray(theta_log10)
    pv[:, _S["gidx"]] = DOSES
    obs = _S["sim"].run(pv, T)
    V = np.array([o["obs_VCAM1"][-1] for o in obs]); P = np.array([o["obs_P"][-1] for o in obs])
    return {("VCAM1", 150): np.log2(V[1] / V[0]), ("VCAM1", 500): np.log2(V[2] / V[0]),
            ("P", 150): np.log2(P[1] / P[0]),     ("P", 500): np.log2(P[2] / P[0])}

def loglik(theta):
    try:
        pred = simulate(theta)
    except Exception:
        return -1e12
    ll = 0.0
    for k, (mu, sd) in DATA.items():
        if not np.isfinite(pred[k]): return -1e12
        ll += -0.5 * ((pred[k] - mu) / sd) ** 2 - np.log(sd * np.sqrt(2 * np.pi))
    return ll

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["M1", "M3"], default="M3")
    ap.add_argument("--niter", type=int, default=3000)
    ap.add_argument("--nchains", type=int, default=3)
    ap.add_argument("--quick", action="store_true", help="skip the convergence check")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    from pydream.core import run_dream
    from pydream.parameters import SampledParam
    from scipy.stats import uniform

    sim, fit_names, fit = setup(a.model)
    print(f"engine: {sim.kind}")
    print(f"{a.model}: fitting {len(fit_names)} parameters: {fit_names}")
    print("data (log2FC, sigma):", {f"{k[0]}@{k[1]}": (round(v[0], 2), round(v[1], 2)) for k, v in DATA.items()})
    t0 = time.time(); ll0 = loglik([np.mean(fit[n]) for n in fit_names]); print(f"one likelihood eval: {time.time()-t0:.2f} s (ll={ll0:.1f})")

    params = [SampledParam(uniform, loc=lo, scale=hi - lo) for lo, hi in (fit[n] for n in fit_names)]
    np.random.seed(a.seed)
    out = os.path.join(HERE, f"dream_{a.model}")
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    sampled, logps = run_dream(params, loglik, niterations=a.niter, nchains=a.nchains, multitry=False,
                               gamma_levels=4, adapt_gamma=False, history_thin=1, model_name=os.path.join(out, "run"), verbose=False)
    print(f"sampling done in {(time.time()-t0)/60:.1f} min")

    burn = a.niter // 2
    samples = np.concatenate([s[burn:] for s in sampled]); lp = np.concatenate([l[burn:] for l in logps]).ravel()
    np.save(os.path.join(out, "samples_log10.npy"), samples); np.save(os.path.join(out, "logp.npy"), lp)

    # ---- summary ----
    best = samples[np.argmax(lp)]
    pred = simulate(best)
    summ = {"model": a.model, "n_params": len(fit_names), "max_loglik": float(lp.max()),
            "AIC": float(2 * len(fit_names) - 2 * lp.max()),
            "best_log10": dict(zip(fit_names, map(float, best))),
            "posterior_quantiles_log10": {n: [float(q) for q in np.quantile(samples[:, i], [0.025, 0.5, 0.975])] for i, n in enumerate(fit_names)},
            "data": {f"{k[0]}@{k[1]}": v for k, v in DATA.items()},
            "best_pred": {f"{k[0]}@{k[1]}": float(v) for k, v in pred.items()}}
    if not a.quick:
        try:
            from pydream.convergence import Gelman_Rubin
            gr = Gelman_Rubin([s[burn:] for s in sampled]); summ["gelman_rubin"] = dict(zip(fit_names, map(float, gr)))
        except Exception as e:
            summ["gelman_rubin"] = f"not computed: {e}"
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)

    print("\nbest fit vs data (log2FC):")
    for k in DATA: print(f"  {k[0]:6s} @ {k[1]:3d}   data {DATA[k][0]:6.2f} ± {DATA[k][1]:.2f}   model {pred[k]:6.2f}")
    print(f"max log-likelihood {lp.max():.2f}   AIC {summ['AIC']:.1f}")
    print("posterior 95% intervals (log10):")
    for n, q in summ["posterior_quantiles_log10"].items():
        width = q[2] - q[0]; tag = "pinned" if width < 1.0 else ("loose" if width < 2.5 else "unidentified")
        print(f"  {n:9s} [{q[0]:6.2f}, {q[2]:6.2f}]  median {q[1]:6.2f}   {tag}")
    print("saved", out)

if __name__ == "__main__":
    main()
