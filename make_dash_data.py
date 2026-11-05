import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bng_setup  # noqa
from pysb.simulator import ScipyOdeSimulator
from pysb.bng import generate_equations
from gev_models import build
import calibrate as C

t0 = time.time()
DOSES = np.arange(0, 601, 20.0)
T = np.linspace(0, 600, 120)
out = {"doses": DOSES.tolist()}

# ---------- structural panel: M1, and M3 over a Km_arr x n_hill grid ----------
def sweep(model, overrides):
    names = [p.name for p in model.parameters]
    base = np.array([p.value for p in model.parameters], float)
    pv = np.tile(base, (len(DOSES), 1))
    for k, v in overrides.items():
        pv[:, names.index(k)] = v
    pv[:, names.index("gEV_0")] = DOSES
    obs = ScipyOdeSimulator(model, tspan=T, compiler="python").run(param_values=pv).observables
    V = np.array([o["obs_VCAM1"][-1] for o in obs]); P = np.array([o["obs_P"][-1] for o in obs])
    return (V / V[0]).tolist(), (P / P[0]).tolist()

M1 = build(with_brake=False); generate_equations(M1)
M3 = build(with_brake=True);  generate_equations(M3)
v, p = sweep(M1, {})
out["m1"] = {"vcam1": v, "p": p}

KMS = [40, 50, 55, 60, 65, 72, 80, 85, 90, 100]
HILLS = [4, 8, 12]
out["kms"] = KMS; out["hills"] = HILLS
grid = {}
for h in HILLS:
    for km in KMS:
        v, p = sweep(M3, {"Km_arr": float(km), "n_hill": float(h)})
        grid[f"{h}_{km}"] = {"vcam1": v, "p": p}
print(f"structural grid done {time.time()-t0:.0f}s")

out["grid"] = grid

# ---------- fit panel ----------
rng = np.random.default_rng(0)
fit = {"data": {}, "models": {}}
for k, (mu, sd) in C.DATA.items():
    fit["data"][f"{k[0]}@{k[1]}"] = {"dose": k[1], "mean": mu, "sd": sd, "obs": k[0]}
for name in ("M1", "M3"):
    C.setup(name)
    S = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"dream_{name}", "summary.json")))
    smp = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"dream_{name}", "samples_log10.npy"))
    lp = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"dream_{name}", "logp.npy"))
    def curves(theta):
        pv = np.tile(C._S["base"], (len(DOSES), 1)); pv[:, C._S["idx"]] = 10 ** theta; pv[:, C._S["gidx"]] = DOSES
        obs = C._S["sim"].run(param_values=pv).observables
        V = np.array([o["obs_VCAM1"][-1] for o in obs]); P = np.array([o["obs_P"][-1] for o in obs])
        return np.log2(V / V[0]).tolist(), np.log2(P / P[0]).tolist()
    bv, bp = curves(smp[np.argmax(lp)])
    draws = []
    for th in smp[rng.choice(len(smp), 25, replace=False)]:
        dv, dp = curves(th); draws.append({"vcam1": dv, "p": dp})
    fit["models"][name] = {
        "best": {"vcam1": bv, "p": bp}, "draws": draws,
        "loglik": S["max_loglik"], "aic": S["AIC"], "n_params": S["n_params"],
        "pred": S["best_pred"], "quantiles": S.get("posterior_quantiles_log10", {}),
        "gelman_rubin": (max(S["gelman_rubin"].values()) if isinstance(S.get("gelman_rubin"), dict) else None),
    }
    print(f"{name} fit curves done {time.time()-t0:.0f}s")
out["fit"] = fit

# round to keep the payload small
def rnd(o):
    if isinstance(o, float): return round(o, 4)
    if isinstance(o, list): return [rnd(x) for x in o]
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    return o
out = rnd(out)
dest = "/tmp/claude-0/-home-claude/7ac629ae-59d3-5434-9a7d-0ef1f22510e9/scratchpad/dash/data.json"
os.makedirs(os.path.dirname(dest), exist_ok=True)
json.dump(out, open(dest, "w"))
print("wrote", dest, os.path.getsize(dest) // 1024, "KB", f"{time.time()-t0:.0f}s")
