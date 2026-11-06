"""
Author-side tool. Regenerates network_M1.py and network_M3.py from gev_models.py.

Running this needs BioNetGen, which needs Perl. Running the MODELS does not: the
generated modules are plain numpy and scipy, which is what lets anyone clone the
repo and run everything without installing Perl.

    python make_network.py

Re-run it whenever gev_models.py changes, then commit the regenerated files.
"""
import os
import sys

import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bng_setup  # noqa: F401
from pysb.bng import generate_equations
from gev_models import build

HEADER = '''"""
Reaction network for {model}, generated from gev_models.py by make_network.py.
Do not edit by hand. Needs only numpy and scipy.

Species are indexed in BioNetGen's order:
{species}
"""
import numpy as np

MODEL = "{model}"
PARAMS = {params}
DEFAULTS = np.array({defaults}, dtype=float)
PARAM_INDEX = {{n: i for i, n in enumerate(PARAMS)}}
# species index -> index of the parameter holding its initial value
INITIALS = {initials}
# species whose concentration is clamped (d/dt forced to zero)
FIXED = {fixed}
# observable name -> list of (species index, coefficient)
OBSERVABLES = {observables}
N_SPECIES = {nspecies}


def rhs(t, y, p):
    return np.array([
{rhs}
    ])


def jac(t, y, p):
    J = np.zeros(({nspecies}, {nspecies}))
{jac}
    return J


def initial_state(p):
    y0 = np.zeros({nspecies})
    for si, pi in INITIALS.items():
        y0[si] = p[pi]
    return y0
'''


def generate(model_name, with_brake):
    m = build(with_brake=with_brake)
    generate_equations(m)
    params = [p.name for p in m.parameters]
    nsp = len(m.species)

    # PySB parameters and expressions are their own Symbol subclasses, so the
    # substitution maps have to be keyed on the objects, not on Symbol(name).
    expr_subs = {e: e.expand_expr(expand_observables=True) for e in m.expressions}

    ysym = [sympy.Symbol("y[%d]" % i) for i in range(nsp)]
    psym = {par: sympy.Symbol("p[%d]" % i) for i, par in enumerate(m.parameters)}
    sspec = [sympy.Symbol("__s%d" % i) for i in range(nsp)]

    def convert(expr):
        e = sympy.sympify(expr)
        for _ in range(6):                       # expressions can nest
            new = e.subs(expr_subs)
            if new == e:
                break
            e = new
        e = e.subs({sspec[i]: ysym[i] for i in range(nsp)})
        e = e.subs(psym)
        return e

    odes = [convert(o) for o in m.odes]
    # the printer must not rename our bracketed symbols
    printer = sympy.printing.pycode

    rhs_lines = []
    for i, o in enumerate(odes):
        rhs_lines.append("        %s,   # d[%s]/dt" % (printer(o), m.species[i]))

    jac_lines = []
    M = sympy.Matrix(odes)
    for i in range(nsp):
        for j in range(nsp):
            d = sympy.simplify(sympy.diff(M[i], sympy.Symbol("y[%d]" % j)))
            if d != 0:
                jac_lines.append("    J[%d, %d] = %s" % (i, j, printer(d)))

    initials = {}
    for ic in m.initials:
        si = None
        for k, s in enumerate(m.species):
            if s.is_equivalent_to(ic.pattern):
                si = k
                break
        initials[si] = params.index(ic.value.name)

    fixed = []
    for ic in m.initials:
        if ic.fixed:
            for k, s in enumerate(m.species):
                if s.is_equivalent_to(ic.pattern):
                    fixed.append(k)

    obs = {o.name: list(zip([int(x) for x in o.species], [int(c) for c in o.coefficients]))
           for o in m.observables}

    src = HEADER.format(
        model=model_name,
        species="\n".join("  %d: %s" % (i, s) for i, s in enumerate(m.species)),
        params=repr(params),
        defaults=repr([float(p.value) for p in m.parameters]),
        initials=repr(initials),
        fixed=repr(sorted(fixed)),
        observables=repr(obs),
        nspecies=nsp,
        rhs="\n".join(rhs_lines),
        jac="\n".join(jac_lines) if jac_lines else "    pass",
    )
    dest = os.path.join(HERE, "network_%s.py" % model_name)
    open(dest, "w", encoding="utf-8").write(src)
    print("wrote %s  (%d species, %d reactions)" % (dest, nsp, len(m.reactions)))


if __name__ == "__main__":
    generate("M1", False)
    generate("M3", True)
