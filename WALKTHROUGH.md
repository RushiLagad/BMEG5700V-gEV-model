# Walkthrough

A guided run of the whole project in 8 to 10 minutes: the structural claim first, then the
calibration. Everything here runs from this folder after `python setup_check.py`, with no
BioNetGen or Perl needed.

## Part 1: the structural claim (3 minutes, no fitting)

Both models have the same 8 species and 12 reactions. M1 is one reversible receptor.
M3 is the same model plus a beta-arrestin2 brake on the occupied receptor. The only
difference between them is `kcat_arr`, which is zero in M1.

```
python demo.py
```

The table prints in about a second, then a figure opens. Read the two "reads as" columns:
M1 is down at both doses, M3 is down and then up. The green curve crosses back over control.

The rate constants here are hand-picked, so this part is about shape, not magnitude. The claim
is structural: M1 cannot produce a sign change across dose at any parameter values, and M3 can.

### Moving the brake threshold

`Km_arr` is the occupancy at which braking reaches half maximum, measured in occupied receptors
out of 100, not in dose units.

```
python demo.py --threshold 60     # the flip moves to a lower dose
python demo.py --threshold 85     # the flip moves right, or off the dose range entirely
```

Set it too low and the 150 dose stops looking protective, which contradicts the data. Set it too
high and the flip never happens. The data therefore bound this parameter from both sides, which
is what makes it estimable.

## Part 2: the calibration (5 minutes)

```
python demo_fit.py --live 150
```

This runs 3 chains of 150 PyDREAM iterations live, about 70 seconds, so the sampling is visible
rather than asserted. A 150-iteration burst is a demonstration, not the fit; the script then
loads the saved 1500-iteration run and reports it.

Four targets are fitted: VCAM1 at 150 and 500, and the flip-protein mean at 150 and 500. The
likelihood is Gaussian on the log2 scale with uniform priors in log10.

### What comes out

| | VCAM1 @150 | VCAM1 @500 | P @150 | P @500 | max log-lik | AIC |
|---|---|---|---|---|---|---|
| Data | -1.24 | -0.48 | -1.79 | +1.62 | | |
| M1 | -0.84 | -0.85 | -0.30 | -0.30 | -27.2 | 66.4 |
| M3 | -0.73 | -0.94 | -1.70 | +1.64 | -10.0 | 40.0 |

Delta AIC is 26 in favour of M3, after paying for four extra parameters. M1's best possible fit
to P is flat at -0.30 for both doses.

### What three doses can and cannot constrain

Measured as shrinkage against the prior, the switch steepness `n_hill` narrows by 65 percent and
the induction amplitude `fold_p` by 56. The rate constants barely move: `kcat_arr` 7 percent,
`kreset` 4. The brake threshold `Km_arr` shows only 28 percent, but that is misleading on its
own, because it correlates -0.90 with the binding affinity `kr`. What the data fix is the dose at
which occupancy crosses the threshold, not either parameter separately.

The line to land: three doses identify the shape of the switch and two parameter combinations,
and say nothing about the individual kinetics. That is a dose-design limit, not a model failure.

### Convergence

Both runs converged over 30000 iterations by 3 chains: Gelman-Rubin 1.003 for M1 and 1.017 for
M3, against a target below 1.1.

`demo_fit.png` is drawn at the end: 40 posterior draws per model, the best fit in bold, data in
black. In the right panel the M1 band is flat through zero while the M3 band crosses it.

## Known misfits

**VCAM1 at 500.** The data say -0.48 and M3 says -0.94. The fit splits the difference between the
two doses rather than reproducing the partial recovery. Either VCAM1 has a recovery mechanism the
model lacks, or the 500 point is noisy; it carries FDR 0.41 in proteomics, and the western blot
establishes only that it is lower, not by how much.

**Hill coefficient.** The posterior prefers n between 6 and 16, steeper than a plain
arrestin-recruitment step should produce. Receptor internalization and recycling, which are not
in the model, are the most likely missing mechanism.

**Placeholder sigma.** The likelihood uses 0.35 log2 units rather than limma's per-protein
standard errors, which have not been exported yet. Every posterior interval depends on it.

## Why only five flip proteins

Tnnc1 is excluded because troponin is absent from smooth muscle, so its 4.54 log2FC at 500 most
likely reflects cardiac tissue in the aortic arch samples; histology is pending. Saa1 and Saa2
are inverted-U but are liver-derived acute-phase proteins rather than artery-wall readouts. The
five that remain, Nceh1, Comp, Serinc3, Aldh1l2 and Sh3bgr, all carry FDR below 1e-3 at both doses.

## Open questions

1. Whether to get the observed steepness from mechanism, by adding receptor internalization and
   recycling, rather than from a large Hill exponent.
2. Whether to model VCAM1's partial recovery at 500 with a second input, or treat that point as
   too noisy to fit given FDR 0.41.
3. Whether `Arr_tot` should be depletable rather than constant. The quasi-steady-state reduction
   assumes arrestin is in excess, but at the top dose a large fraction of receptors is braked.
