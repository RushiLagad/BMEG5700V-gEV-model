# Calibration of M1 vs M3 against the real data

Data: log2 fold change vs control at 150 and 500 mg/kg/day, from limma (Combined report, Aug 20).
  VCAM1:        -1.24 @150,  -0.48 @500
  flip protein: mean of 5 anchors (Nceh1, Comp, Serinc3, Aldh1l2, Sh3bgr): -1.79 @150, +1.62 @500
  Tnnc1 excluded pending histology. Saa1/Saa2 excluded (acute-phase, liver-derived).
  Error model: Gaussian on log2 scale, sigma = max(spread across anchors, 0.35). PLACEHOLDER until
  limma standard errors are exported.

PyDREAM: 3 chains x 30000 iterations per model, second half kept (45000 samples each).
Uniform priors in log10, ranges in calibrate.py. Both runs converged:
Gelman-Rubin max 1.003 for M1 and 1.017 for M3, against a target below 1.1.

| | parameters | max log-lik | AIC | P@150 | P@500 | VCAM1@150 | VCAM1@500 |
|---|---|---|---|---|---|---|---|
| data | | | | -1.79 | +1.62 | -1.24 | -0.48 |
| M1 one receptor | 6 | -27.2 | 66.4 | -0.30 | -0.30 | -0.84 | -0.85 |
| M3 receptor + brake | 10 | -10.0 | 40.0 | -1.70 | +1.64 | -0.73 | -0.94 |

Delta AIC = 26 in favour of M3, even after paying for 4 extra parameters.
M1's best compromise is a flat P: it cannot go down at 150 and up at 500 at any parameter values.
That is a structural statement about one saturating receptor, not a statement about this fit.

## What the three doses actually constrain (M3)

Measured as shrinkage: how much narrower the 95% posterior is than the prior it started from.

| parameter | 95% posterior (log10) | shrinkage | reading |
|---|---|---|---|
| n_hill | [0.78, 1.20] | 65% | constrained. Hill coefficient 6 to 16 |
| fold_p | [0.83, 1.49] | 56% | constrained. 7x to 31x max induction of P by the brake arm |
| kr | [1.95, 3.49] | 49% | constrained |
| Km_arr | [0.92, 1.99] | 28% | weak on its own, see the trade-off below |
| KI_v | [0.15, 2.47] | 23% | weak |
| KA_p | [-0.46, 1.98] | 19% | essentially the prior |
| KI_p | [-0.45, 2.04] | 17% | essentially the prior |
| ksig | [-0.97, 1.30] | 9% | essentially the prior |
| kcat_arr | [-1.83, 1.44] | 7% | essentially the prior |
| kreset | [-1.45, 0.95] | 4% | essentially the prior |

Parameters that trade off against each other, so only the combination is identified:

| pair | correlation |
|---|---|
| KI_v vs KI_p | +0.90 |
| kr vs Km_arr | -0.90 |
| ksig vs KI_p | +0.83 |
| ksig vs KI_v | +0.82 |

Reading: three doses identify the SHAPE of the switch (how steep it is, how much induction it
produces) and two parameter combinations. They say almost nothing about the individual rate
constants. The brake threshold `Km_arr` looks weakly constrained on its own only because it
trades off almost perfectly against the binding affinity `kr`: what the data fix is the dose at
which occupancy crosses the threshold, not either number separately. That is the honest headline,
and it is an argument for more doses rather than for a different model.

Note this corrects an earlier, unconverged run which reported `Km_arr` as pinned. It was not;
that was a sampling artifact.

## Caveats, in order of importance
1. The Hill coefficient wants to be high (6 to 16). Open question: whether that is defensible for
   arrestin recruitment or whether a slower process (receptor internalisation / recycling) should
   carry part of the switch. This is the main modeling question the fit raises.
2. VCAM1 is fitted as a compromise (-0.73 / -0.94 vs -1.24 / -0.48). The data say VCAM1 partly
   RECOVERS at the high dose; neither model does that. Possible fix: let the braked receptor
   contribute to VCAM1 recovery, or treat the 500 value as the western blot "down" only, given
   that it carries FDR 0.41 in proteomics.
3. sigma is a placeholder. Export limma's per-protein SE (sqrt(s2.post) * stdev.unscaled) and
   replace SIGMA_FLOOR with per-point sigmas. Every interval above depends on this choice.
4. Dose units: the receptor Kd is fitted in nominal mg/kg/day units with no plasma PK step. Fine
   for shape; not interpretable as a binding constant.
5. The five flip-protein anchors are averaged into one readout. If histology shows cardiac
   contamination in the aortic arch samples, that readout changes and so does the fit.

## Files
  targets.csv                     the 14 numbers and what is used
  calibrate.py                    the PyDREAM harness (python calibrate.py --model M3 --niter 30000)
  dream_M1/, dream_M3/            samples_log10.npy, logp.npy, summary.json, plus _full reference copies
  demo_fit.py                     prints this table and draws demo_fit.png
  demo_fit.png                    posterior predictive against the data
