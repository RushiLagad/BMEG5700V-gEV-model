"""
gamma-EV dose-response models (BMEG 5700V project), v2, revised after code review.

Two candidate wirings of the direct (receptor) route, written in PySB.
Units are deliberately arbitrary: the goal at this stage is shape, not numbers.

Changes from v1 (his two comments):
  1. Binding is a proper reversible rule,  gEV + CaSR | gEV%CaSR  (kf, kr).
     v1 wrote  gEV + CaSR >> gEV + CaSR%gEV , which is a synthesis rule, not binding.
     The dose is held constant with a FIXED initial condition instead.
  2. The brake is a quasi-steady-state (Michaelis-Menten / Hill) reduction:
     arrestin is the "enzyme" in excess (Arr_tot), the OCCUPIED RECEPTOR is the
     substrate, and the Hill term sits on that substrate:
         flux = kcat * Arr_tot * [occ]^n / (Km^n + [occ]^n)
     v1 multiplied a Hill term in free gEV by mass action in two other species.
     Because PySB multiplies a unimolecular rule's rate by its reactant, the rate
     expression on the rule is flux / [occ].

M1  one saturating receptor, no brake
    gEV + CaSR <-> gEV%CaSR             reversible binding
    gEV%CaSR   -> gEV%CaSR + Sig        active receptor makes a calming signal
    Sig        -> 0
    Sig represses VCAM1 and P synthesis

M3  M1 plus a beta-arrestin2 brake
    gEV%CaSR(arr=off) -> gEV%CaSR(arr=on)   Hill-type, see above (the brake)
    gEV%CaSR(arr=on)  -> gEV + CaSR         reset / recycling; releases the ligand
    braked receptor makes no Sig
    braked receptor makes S2, which induces P
"""
from pysb import Model, Monomer, Parameter, Rule, Initial, Observable, Expression


def build(with_brake: bool):
    m = Model(name="M3_receptor_plus_brake" if with_brake else "M1_one_receptor", _export=False)

    def add(c):
        m.add_component(c)
        return c

    # ---- species ---------------------------------------------------------
    gEV   = add(Monomer("gEV", ["b"], _export=False))
    CaSR  = add(Monomer("CaSR", ["b", "arr"], {"arr": ["off", "on"]}, _export=False))
    Sig   = add(Monomer("Sig", _export=False))
    S2    = add(Monomer("S2", _export=False))
    VCAM1 = add(Monomer("VCAM1", _export=False))
    P     = add(Monomer("P", _export=False))

    # ---- parameters (arbitrary units) --------------------------------------
    prm = {}
    for name, val in [
        ("gEV_0", 0.0),                       # dose; held fixed (see Initial below)
        ("CaSR_0", 100.0),
        ("kf", 1.0), ("kr", 150.0),           # Kd = kr/kf = 150 dose units
        ("ksig", 1.0), ("kdsig", 1.0),
        ("kv_basal", 10.0), ("kv_deg", 1.0), ("KI_v", 4.0),
        ("kp_basal", 10.0), ("kp_deg", 1.0), ("KI_p", 30.0),
        # brake (Michaelis-Menten / Hill reduction); kcat_arr = 0 gives M1
        ("kcat_arr", 0.3 if with_brake else 0.0),
        ("Arr_tot", 100.0),
        ("Km_arr", 72.0),                     # occupied-receptor level at half-max braking
        ("n_hill", 12.0),
        ("kreset", 0.5),
        ("ks2", 1.0), ("kds2", 1.0), ("KA_p", 10.0), ("fold_p", 8.0),
        ("VCAM1_0", 10.0), ("P_0", 10.0),
    ]:
        prm[name] = add(Parameter(name, val, _export=False))

    # ---- initial conditions ----------------------------------------------
    m.add_initial(Initial(gEV(b=None), prm["gEV_0"], fixed=True, _export=False))   # clamped dose
    m.add_initial(Initial(CaSR(b=None, arr="off"), prm["CaSR_0"], _export=False))
    m.add_initial(Initial(VCAM1(), prm["VCAM1_0"], _export=False))
    m.add_initial(Initial(P(), prm["P_0"], _export=False))

    # ---- observables -----------------------------------------------------
    occ    = add(Observable("obs_occ", CaSR(b=1, arr="off") % gEV(b=1), _export=False))
    braked = add(Observable("obs_braked", CaSR(b=1, arr="on") % gEV(b=1), _export=False))
    sig    = add(Observable("obs_Sig", Sig(), _export=False))
    s2     = add(Observable("obs_S2", S2(), _export=False))
    add(Observable("obs_free_R", CaSR(b=None), _export=False))
    add(Observable("obs_VCAM1", VCAM1(), _export=False))
    add(Observable("obs_P", P(), _export=False))

    # ---- 1. reversible binding (Harris comment 1) --------------------------
    add(Rule("bind", gEV(b=None) + CaSR(b=None, arr="off") | gEV(b=1) % CaSR(b=1, arr="off"),
             prm["kf"], prm["kr"], _export=False))

    # active, un-braked receptor produces the calming signal
    add(Rule("make_sig", CaSR(b=1, arr="off") % gEV(b=1) >> CaSR(b=1, arr="off") % gEV(b=1) + Sig(),
             prm["ksig"], _export=False))
    add(Rule("sig_decay", Sig() >> None, prm["kdsig"], _export=False))

    # ---- 2. the brake as a Hill-type QSSA reduction (Harris comment 2) -----
    # flux = kcat * Arr_tot * occ^n / (Km^n + occ^n); PySB multiplies by [occ], so divide once.
    n = prm["n_hill"]
    brake_rate = add(Expression("brake_rate",
                                prm["kcat_arr"] * prm["Arr_tot"] * occ ** (n - 1) / (prm["Km_arr"] ** n + occ ** n),
                                _export=False))
    add(Rule("brake_on", CaSR(b=1, arr="off") % gEV(b=1) >> CaSR(b=1, arr="on") % gEV(b=1),
             brake_rate, _export=False))
    # reset: braked receptor recycles and releases its ligand (his rule 3)
    add(Rule("brake_reset", CaSR(b=1, arr="on") % gEV(b=1) >> CaSR(b=None, arr="off") + gEV(b=None),
             prm["kreset"], _export=False))
    add(Rule("make_s2", CaSR(b=1, arr="on") % gEV(b=1) >> CaSR(b=1, arr="on") % gEV(b=1) + S2(),
             prm["ks2"], _export=False))
    add(Rule("s2_decay", S2() >> None, prm["kds2"], _export=False))

    # ---- readouts --------------------------------------------------------
    v_syn = add(Expression("v_syn", prm["kv_basal"] * prm["KI_v"] / (prm["KI_v"] + sig), _export=False))
    add(Rule("vcam_syn", None >> VCAM1(), v_syn, _export=False))
    add(Rule("vcam_deg", VCAM1() >> None, prm["kv_deg"], _export=False))

    p_syn = add(Expression("p_syn", prm["kp_basal"] * (prm["KI_p"] / (prm["KI_p"] + sig))
                           * (1 + prm["fold_p"] * s2 / (prm["KA_p"] + s2)), _export=False))
    add(Rule("p_syn_rule", None >> P(), p_syn, _export=False))
    add(Rule("p_deg", P() >> None, prm["kp_deg"], _export=False))
    return m


M1 = build(with_brake=False)
M3 = build(with_brake=True)
