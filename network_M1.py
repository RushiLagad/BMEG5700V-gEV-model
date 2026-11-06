"""
Reaction network for M1, generated from gev_models.py by make_network.py.
Do not edit by hand. Needs only numpy and scipy.

There is deliberately no analytic Jacobian here. Differentiating the Hill term leaves
x**n/x, which is 0/0 at zero occupancy, and the tidy symbolic forms put x back under
the fraction. A finite-difference Jacobian costs about 3 percent more time and cannot
produce a NaN, so the solver computes its own.

Species are indexed in BioNetGen's order:
  0: gEV(b=None)
  1: CaSR(b=None, arr='off')
  2: VCAM1()
  3: P()
  4: CaSR(b=1, arr='off') % gEV(b=1)
  5: Sig()
  6: CaSR(b=1, arr='on') % gEV(b=1)
  7: S2()
"""
import numpy as np

MODEL = "M1"
PARAMS = ['gEV_0', 'CaSR_0', 'kf', 'kr', 'ksig', 'kdsig', 'kv_basal', 'kv_deg', 'KI_v', 'kp_basal', 'kp_deg', 'KI_p', 'kcat_arr', 'Arr_tot', 'Km_arr', 'n_hill', 'kreset', 'ks2', 'kds2', 'KA_p', 'fold_p', 'VCAM1_0', 'P_0']
DEFAULTS = np.array([0.0, 100.0, 1.0, 150.0, 1.0, 1.0, 10.0, 1.0, 4.0, 10.0, 1.0, 30.0, 0.0, 100.0, 72.0, 12.0, 0.5, 1.0, 1.0, 10.0, 8.0, 10.0, 10.0], dtype=float)
PARAM_INDEX = {n: i for i, n in enumerate(PARAMS)}
# species index -> index of the parameter holding its initial value
INITIALS = {0: 0, 1: 1, 2: 21, 3: 22}
# species whose concentration is clamped (d/dt forced to zero)
FIXED = [0]
# observable name -> list of (species index, coefficient)
OBSERVABLES = {'obs_occ': [(4, 1)], 'obs_braked': [(6, 1)], 'obs_Sig': [(5, 1)], 'obs_S2': [(7, 1)], 'obs_free_R': [(1, 1)], 'obs_VCAM1': [(2, 1)], 'obs_P': [(3, 1)]}
N_SPECIES = 8


def rhs(t, y, p):
    return np.array([
        0,   # d[gEV(b=None)]/dt
        p[16]*y[6] - p[2]*y[0]*y[1] + p[3]*y[4],   # d[CaSR(b=None, arr='off')]/dt
        p[6]*p[8]/(p[8] + y[5]) - p[7]*y[2],   # d[VCAM1()]/dt
        -p[10]*y[3] + p[11]*p[9]*(p[20]*y[7]/(p[19] + y[7]) + 1)/(p[11] + y[5]),   # d[P()]/dt
        -p[12]*p[13]*y[4]*y[4]**(p[15] - 1)/(p[14]**p[15] + y[4]**p[15]) + p[2]*y[0]*y[1] - p[3]*y[4],   # d[CaSR(b=1, arr='off') % gEV(b=1)]/dt
        p[4]*y[4] - p[5]*y[5],   # d[Sig()]/dt
        p[12]*p[13]*y[4]*y[4]**(p[15] - 1)/(p[14]**p[15] + y[4]**p[15]) - p[16]*y[6],   # d[CaSR(b=1, arr='on') % gEV(b=1)]/dt
        p[17]*y[6] - p[18]*y[7],   # d[S2()]/dt
    ])


def initial_state(p):
    y0 = np.zeros(8)
    for si, pi in INITIALS.items():
        y0[si] = p[pi]
    return y0
