# γ-EV dose-response models (PySB)

Mechanistic models of how γ-glutamyl valine (γ-EV, a dipeptide from Great Northern beans)
produces a non-monotonic protein response in the aorta of ApoE-deficient mice.
Semester project for BMEG 4700V/5700V, Introduction to Computational and Systems Biology,
University of Arkansas, Fall 2026. Rushikesh R. Lagad, Goswami Lab.

**The question.** At 150 mg/kg/day the aortic proteome moves one way; at 500 mg/kg/day a set of
proteins moves back the other way, while VCAM1 stays down at both doses. One saturating receptor
cannot do that. The models ask whether a β-arrestin2 brake on the occupied receptor can, and the
calibration asks how much of it the three available doses can actually constrain.

**Data provenance.** `targets.csv` holds limma log2 fold changes from an aortic LFQ proteomics
experiment run in collaboration with Dr. Kaustav Majumder's lab (University of Nebraska-Lincoln).
These results are unpublished. The dose labels are nominal in the model, which fits dose in
arbitrary units; see the caveats below before reusing any number here.

**Interactive explorer.** An interactive version of the two panels below, with a slider for the
brake threshold and the fitted posteriors drawn against the data:
https://claude.ai/artifact/GeYWW9LWJCcEsR8hXijkyM (viewable by people it has been shared with).

Two candidate wirings of the direct receptor route, defined in PySB and simulated from the
pre-generated networks. Both models have 8 species, 12 reactions, 8 ODEs.

| File | What it is |
|---|---|
| `gev_models.py` | The two models. `build(with_brake=False)` is M1, `build(with_brake=True)` is M3. Same code; M1 just has `kcat_arr` set to zero. |
| `run_dose_response.py` | Runs each model to steady state across 12 doses and plots VCAM1, the flip protein P, and the receptor states. |
| `demo.py` | Fast structural demo. `--threshold N` moves the brake threshold `Km_arr`, `--hill N` the Hill coefficient, `--no-plot` to skip the figure. |
| `targets.csv` | The real limma log2 fold changes used as fit targets, with FDR, source, and why each protein is in or out. |
| `calibrate.py` | PyDREAM calibration of either model against `targets.csv`. |
| `demo_fit.py` | Live demo of the calibration. `--live N` runs a short visible PyDREAM burst first, then loads the saved full run. |
| `CHANGES_after_review.md` | Each code-review comment mapped to what changed in the model. |
| `FIT_RESULTS.md` | The fitted numbers, what is identifiable, and the caveats. |
| `WALKTHROUGH.md` | A guided run of the whole thing in 8 to 10 minutes. |
| `engine.py` | One simulation interface. Uses the pre-generated networks by default, PySB if you ask for it. |
| `network_M1.py`, `network_M3.py` | The reaction networks generated from the PySB models. Plain numpy and scipy, no BioNetGen. |
| `setup_check.py` | One-time install, the PyDREAM NumPy 2 patch, and a smoke test. |
| `requirements.txt` | What you need to run it. |
| `dream_M1/`, `dream_M3/` | Saved posterior samples, log-likelihoods, and fit summaries. |
| `make_network.py` | Author-side. Regenerates the network files from `gev_models.py`. Needs BioNetGen. |
| `verify_network.py` | Author-side. Checks the generated networks still match the PySB models. |
| `bng_setup.py`, `fix_bng.py` | Author-side helpers for getting BioNetGen working on Windows. |

## Setup (one time)
```
git clone https://github.com/RushiLagad/BMEG5700V-gEV-model.git
cd BMEG5700V-gEV-model
python setup_check.py
```
That installs numpy, scipy, pandas, matplotlib and pydream if they are missing, applies the
PyDREAM NumPy 2 patch, and smoke-tests both models and the saved fits. Python 3.10 or newer.

**You do not need BioNetGen or Perl.** `gev_models.py` is the PySB source of truth for both
models, but its reaction network is pre-generated into `network_M1.py` and `network_M3.py`,
which contain only numpy and scipy. That is what everything runs on, and it is also about ten
times faster than going through PySB.

To work from the PySB models directly instead, install `pysb` and `bionetgen` (and Perl, on
Windows) and either set `GEV_ENGINE=pysb` or pass `engine="pysb"` to `engine.load`. After
editing `gev_models.py`, run `python make_network.py` to regenerate, then `python
verify_network.py` to confirm the two agree. They currently match to 2.4e-8.

## Run
```
python demo.py                          # structural claim, about a second
python run_dose_response.py             # full dose sweep figure
python demo_fit.py                      # load the saved fit, print M1 vs M3, draw demo_fit.png
python demo_fit.py --live 150           # 70 s of visible sampling first
python calibrate.py --model M3 --niter 30000 --nchains 3   # the run behind the numbers below
```

## What the fit says (real targets, PyDREAM, 3 chains x 30000, both converged)
| | VCAM1 @150 | VCAM1 @500 | P @150 | P @500 | max log-lik | AIC |
|---|---|---|---|---|---|---|
| Data (log2FC) | -1.24 | -0.48 | -1.79 | +1.62 | | |
| M1 one receptor | -0.84 | -0.85 | -0.30 | -0.30 | -27.2 | 66.4 |
| M3 receptor + brake | -0.73 | -0.94 | -1.70 | +1.64 | -10.0 | 40.0 |

Delta AIC = 26 in favour of M3, after paying for its 4 extra parameters. M1's best possible fit
to P is flat: one saturating receptor cannot change sign across dose at any parameter values.

Gelman-Rubin is 1.003 for M1 and 1.017 for M3, so both runs have converged.

Of M3's 10 fitted parameters, three doses constrain the shape of the switch and little else.
Measured as shrinkage against the prior: `n_hill` 65%, `fold_p` 56%, `kr` 49%, `Km_arr` 28%, and
the remaining six between 4% and 23%. `Km_arr` correlates -0.90 with `kr`, so what the data fix
is the dose at which occupancy crosses the brake threshold, not either parameter alone. Full
table and the other trade-offs in `FIT_RESULTS.md`.

## Honest caveats
* The brake is keyed to the occupied receptor through a Michaelis-Menten (Hill) reduction with
  beta-arrestin2 as the enzyme held at constant total, following code review. An earlier version
  keyed it to the clamped peptide level, which was wrong; see `CHANGES_after_review.md`.
* Dose enters as `gEV_0` in the same arbitrary units as `kf`/`kr`, so the receptor Kd is fitted
  rather than assumed. The mg/kg/day labels on the x axis are nominal.
* The likelihood sigma is a 0.35 log2 placeholder, not limma's per-protein standard errors.
  Every posterior interval depends on that choice.
* The posterior wants a Hill coefficient of 6 to 16, steeper than a plain recruitment step should
  give. Receptor internalization and recycling are the most likely missing mechanism.
* VCAM1 recovers partway at 500 in the data (-0.48) and the model does not reproduce it (-0.94).
* The individual rate constants are not identified by three doses. The model ranking is robust;
  the kinetics are not, and should not be quoted as estimates.

## Next
1. Export limma per-protein SEs and replace `SIGMA_FLOOR`.
2. Histology on the aortic arch samples to settle whether the flip proteins are cardiac contamination.
3. Add the gut/acetate branch (M2, M4) as a structural "can it fit at all" test.
4. More doses, if a follow-up study is possible: the trade-offs above are a dose-design problem.

## Licence
MIT, see `LICENSE`. The code is free to reuse. The unpublished fold changes in `targets.csv`
are the collaboration's data, not covered by that licence; ask before building on them.
