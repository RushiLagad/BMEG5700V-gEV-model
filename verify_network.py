"""
Check that the generated networks still match the PySB models they came from.

Needs BioNetGen (and so Perl), which is why it is not part of setup_check.py. Run it
after changing gev_models.py and regenerating with make_network.py.

    python verify_network.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from engine import load

T = np.linspace(0, 600, 60)
DOSES = [0.0, 150.0, 500.0]
rng = np.random.default_rng(0)
worst = 0.0
fail = False

for name in ("M1", "M3"):
    net = load(name, engine="network")
    pysb = load(name, engine="pysb")
    assert net.names == pysb.names, "parameter order differs between the two engines"

    # defaults, then random perturbations of the fitted parameters
    trials = [net.defaults.copy()]
    vary = [n for n in ("kr", "ksig", "KI_v", "KI_p", "fold_p", "KA_p",
                        "kcat_arr", "Km_arr", "n_hill", "kreset") if n in net.index]
    for _ in range(6):
        p = net.defaults.copy()
        for v in vary:
            p[net.index[v]] *= 10 ** rng.uniform(-1, 1)
        trials.append(p)

    for ti, base in enumerate(trials):
        pv = np.tile(base, (len(DOSES), 1))
        pv[:, net.index["gEV_0"]] = DOSES
        a = net.run(pv, T)
        b = pysb.run(pv, T)
        for i in range(len(DOSES)):
            for obs in ("obs_VCAM1", "obs_P", "obs_occ", "obs_braked"):
                x, y = a[i][obs][-1], b[i][obs][-1]
                denom = max(abs(x), abs(y), 1e-12)
                rel = abs(x - y) / denom
                worst = max(worst, rel)
                if rel > 1e-6:
                    fail = True
                    print("MISMATCH %s trial %d dose %g %s: network %.10g vs pysb %.10g"
                          % (name, ti, DOSES[i], obs, x, y))
    print("%s checked over %d parameter sets x %d doses" % (name, len(trials), len(DOSES)))

print("\nworst relative difference: %.2e" % worst)
print("FAILED" if fail else "OK, the generated networks match the PySB models")
sys.exit(1 if fail else 0)
