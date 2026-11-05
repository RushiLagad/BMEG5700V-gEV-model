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

Two candidate wirings of the direct receptor route, built in PySB and simulated through
BioNetGen. Both models have 8 species, 12 reactions, 8 ODEs.

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
| `setup_check.py` | One-time install, the PyDREAM NumPy 2 patch, and a smoke test. |
| `requirements.txt` | The exact package versions this was run and calibrated with. |
| `bng_setup.py` | Sets `BNGPATH` so BioNetGen finds its own Perl modules on Windows. |
| `fix_bng.py` | Diagnostic for the BioNetGen / perl setup if `demo.py` fails. |
| `dream_M1/`, `dream_M3/` | Saved posterior samples, log-likelihoods, and fit summaries. |

## Setup (one time)
```
python setup_check.py
```
That installs `pysb`, `bionetgen`, `pydream`, `pandas`, `matplotlib` and `scipy` if missing,
applies the PyDREAM NumPy 2 patch, and runs a smoke test.

Doing it by hand instead: `pip install pysb bionetgen matplotlib pandas pydream`, then either
`pip install "numpy<2"` or edit `pydream/Dream.py` lines 468 and 523, changing `[0])` to `[0][0])`
(PyDREAM 2.0.0 calls `int()` on a one-element array, which NumPy 2 refuses).

On Windows, if BioNetGen dies with `Can't locate BNGUtils.pm in @INC`, the perl first on PATH is
the MSYS / Git one, whose `FindBin` cannot resolve a `C:\` script path. `bng_setup.py` fixes this
by setting `BNGPATH`, and the entry scripts import it. `python fix_bng.py` diagnoses the rest.

## Run
```
python demo.py                          # structural claim, about a second
python run_dose_response.py             # full dose sweep figure
python demo_fit.py                      # load the saved fit, print M1 vs M3, draw demo_fit.png
python demo_fit.py --live 150           # 70 s of visible sampling first
python calibrate.py --model M3 --niter 5000 --nchains 3    # the real run, tens of minutes
```

## What the fit says (real targets, PyDREAM, 3 chains x 1500)
| | VCAM1 @150 | VCAM1 @500 | P @150 | P @500 | max log-lik | AIC |
|---|---|---|---|---|---|---|
| Data (log2FC) | -1.24 | -0.48 | -1.79 | +1.62 | | |
| M1 one receptor | -0.87 | -0.89 | -0.31 | -0.31 | -27.2 | 66.4 |
| M3 receptor + brake | -0.70 | -0.91 | -1.71 | +1.43 | -10.2 | 40.4 |

Delta AIC = 26 in favour of M3, after paying for its 4 extra parameters. M1's best possible fit
to P is flat: one saturating receptor cannot change sign across dose at any parameter values.

Of M3's 10 fitted parameters, three doses pin down three: `fold_p`, `Km_arr`, `n_hill`.
`KI_p`, `KA_p` and `kcat_arr` are unidentified from this data.

## Honest caveats
* The brake is keyed to the occupied receptor through a Michaelis-Menten (Hill) reduction with
  beta-arrestin2 as the enzyme held at constant total, following code review. An earlier version
  keyed it to the clamped peptide level, which was wrong; see `CHANGES_after_review.md`.
* Dose enters as `gEV_0` in the same arbitrary units as `kf`/`kr`, so the receptor Kd is fitted
  rather than assumed. The mg/kg/day labels on the x axis are nominal.
* The M3 run has not converged: Gelman-Rubin max 1.61 against a target of 1.1. The AIC gap is
  large enough that the model ranking stands, but no posterior interval should be quoted yet.
* The likelihood sigma is a 0.35 log2 placeholder, not limma's per-protein standard errors.
* The posterior wants a Hill coefficient of 7 to 15, steeper than a plain recruitment step should
  give. Probably a missing mechanism, most likely receptor internalization and recycling.
* VCAM1 recovers partway at 500 in the data (-0.48) and the model does not reproduce it (-0.91).

## Next
1. Run `calibrate.py --model M3 --niter 5000` to convergence before quoting intervals.
2. Export limma per-protein SEs and replace `SIGMA_FLOOR`.
3. Histology on the aortic arch samples to settle whether the flip proteins are cardiac contamination.
4. Add the gut/acetate branch (M2, M4) as a structural "can it fit at all" test.
