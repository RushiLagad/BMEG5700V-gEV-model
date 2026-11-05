# First calibration of M1 vs M3 against the real data (Nov 5)

Data: log2 fold change vs control at 150 and 500 mg/kg/day, from limma (Combined report, Aug 20).
  VCAM1:        -1.24 @150,  -0.48 @500
  flip protein: mean of 5 anchors (Nceh1, Comp, Serinc3, Aldh1l2, Sh3bgr): -1.79 @150, +1.62 @500
  Tnnc1 excluded pending histology. Saa1/Saa2 excluded (acute-phase, liver-derived).
  Error model: Gaussian on log2 scale, sigma = max(spread across anchors, 0.35). PLACEHOLDER until
  limma standard errors are exported.

PyDREAM: 3 chains x 1500 iterations, second half kept. Uniform priors in log10 (ranges in calibrate.py).

| | parameters | max log-lik | AIC | P@150 | P@500 | VCAM1@150 | VCAM1@500 |
|---|---|---|---|---|---|---|---|
| data | | | | -1.79 | +1.62 | -1.24 | -0.48 |
| M1 one receptor | 6 | -27.2 | 66.4 | -0.31 | -0.31 | -0.87 | -0.89 |
| M3 receptor + brake | 10 | -10.2 | 40.4 | -1.71 | +1.43 | -0.70 | -0.91 |

Delta AIC = 26 in favour of M3, even after paying for 4 extra parameters.
M1's best compromise is a flat P: it cannot go down at 150 and up at 500 at any parameter values.

## Which parameters the data pins down (M3, 95% posterior interval, log10)
  pinned:        fold_p [0.80, 1.49]  (6x to 31x max induction of P by the brake arm)
                 Km_arr [1.07, 1.98]  (brake half-max at 12 to 95 occupied receptors of 100)
                 n_hill [0.82, 1.19]  (Hill coefficient 7 to 15: steep)
  loose:         kr (Kd), ksig, KI_v, kreset
  unidentified:  KI_p, KA_p, kcat_arr

Reading: three doses pin the SHAPE parameters (how strong the induction is, where the brake sits,
how sharp it is) and say almost nothing about the rate constants. That is the honest headline.

## Caveats, in order of importance
1. Gelman-Rubin for M3 is 1.2 to 1.6 on several parameters: NOT converged at 1500 iterations.
   Run 5000+ iterations (about 30 min on a laptop) before quoting intervals. M1 is converged (GR < 1.1).
2. The Hill coefficient wants to be high (7 to 15). Open question: whether that is defensible for
   arrestin recruitment or whether a slower process (receptor internalisation / recycling) should
   carry part of the switch. This is the main modeling question the fit raises.
3. VCAM1 is fitted as a compromise (-0.70 / -0.91 vs -1.24 / -0.48). The data says VCAM1 partly
   RECOVERS at the high dose; neither model does that well. Possible fix: let the braked receptor
   contribute a little to VCAM1 recovery, or treat the 500 value as the western blot "down" only.
4. sigma is a placeholder. Export limma's per-protein SE (sqrt(s2.post) * stdev.unscaled) and
   replace SIGMA_FLOOR with per-point sigmas.
5. Dose units: the receptor Kd is fitted in mg/kg/day units (no plasma PK step). Fine for shape;
   not interpretable as a binding constant.

## Files
  targets.csv                     the 14 numbers and what is used
  calibrate.py                    the PyDREAM harness (python calibrate.py --model M3 --niter 5000)
  dream_M1/, dream_M3/            samples_log10.npy, logp.npy, summary.json
  posterior_predictive_M1_vs_M3.png
