"""
One simulation interface for both models, with no BioNetGen or Perl needed at run time.

By default this integrates the pre-generated reaction networks in network_M1.py and
network_M3.py, which contain nothing but numpy and scipy. That is what lets anyone
clone the repo and run the models.

If you have BioNetGen working and want to go through PySB instead, either pass
engine="pysb" or set the environment variable GEV_ENGINE=pysb. The two agree to
about 1e-9; verify_network.py checks that.

    from engine import load
    sim = load("M3")
    pv = np.tile(sim.defaults, (3, 1))
    pv[:, sim.index["gEV_0"]] = [0, 150, 500]
    out = sim.run(pv, np.linspace(0, 600, 60))
    out[1]["obs_P"][-1]
"""
import os
import sys

import numpy as np
from scipy.integrate import odeint

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


class _NetworkSim:
    """Integrates a generated network module."""

    kind = "network"

    def __init__(self, model_name):
        mod = __import__("network_%s" % model_name)
        self.mod = mod
        self.model_name = model_name
        self.names = list(mod.PARAMS)
        self.defaults = np.array(mod.DEFAULTS, dtype=float)
        self.index = dict(mod.PARAM_INDEX)
        self.observables = list(mod.OBSERVABLES)

    def run(self, param_values, tspan):
        pv = np.atleast_2d(np.asarray(param_values, dtype=float))
        tspan = np.asarray(tspan, dtype=float)
        mod = self.mod
        out = []
        for p in pv:
            y0 = mod.initial_state(p)
            y = odeint(lambda yy, tt: mod.rhs(tt, yy, p), y0, tspan,
                       Dfun=lambda yy, tt: mod.jac(tt, yy, p),
                       rtol=1e-8, atol=1e-10, mxstep=20000)
            rec = {}
            for name, terms in mod.OBSERVABLES.items():
                acc = np.zeros(len(tspan))
                for si, coeff in terms:
                    acc += coeff * y[:, si]
                rec[name] = acc
            out.append(rec)
        return out


class _PysbSim:
    """Goes through PySB and BioNetGen. Needs Perl on Windows."""

    kind = "pysb"

    def __init__(self, model_name):
        import bng_setup  # noqa: F401
        from pysb.bng import generate_equations
        from pysb.simulator import ScipyOdeSimulator
        from gev_models import build
        self.model_name = model_name
        m = build(with_brake=(model_name == "M3"))
        generate_equations(m)
        self.model = m
        self.names = [p.name for p in m.parameters]
        self.defaults = np.array([p.value for p in m.parameters], dtype=float)
        self.index = {n: i for i, n in enumerate(self.names)}
        self.observables = [o.name for o in m.observables]
        self._ScipyOdeSimulator = ScipyOdeSimulator

    def run(self, param_values, tspan):
        # tight tolerances so this agrees with the generated network; PySB's
        # defaults are looser and drift by about 5e-5 on slow-converging runs
        sim = self._ScipyOdeSimulator(self.model, tspan=np.asarray(tspan, float),
                                      compiler="python",
                                      integrator_options={"rtol": 1e-10, "atol": 1e-12})
        return sim.run(param_values=np.atleast_2d(param_values)).observables


def load(model_name, engine=None):
    """model_name is "M1" or "M3". engine is "network" (default), "pysb", or None."""
    choice = engine or os.environ.get("GEV_ENGINE", "network")
    if choice == "pysb":
        return _PysbSim(model_name)
    try:
        return _NetworkSim(model_name)
    except ImportError:
        print("network_%s.py is missing, falling back to PySB and BioNetGen. "
              "Run make_network.py to regenerate it." % model_name)
        return _PysbSim(model_name)
