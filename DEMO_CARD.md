# Demo card: M1 vs M3, structure then fit (8 to 10 minutes)

v3, after Dr. Harris's review. Two halves: the structural claim (fast, no fitting) and the
calibration (PyDREAM against the real log2 fold changes).

## Before you walk in
```
pip install pysb bionetgen matplotlib pandas pydream
cd gev_model
python demo.py --no-plot        # prints the table in about a second
python demo_fit.py --no-plot    # prints the M1 vs M3 comparison in about two seconds
```
If both print, you are ready. If PyDREAM errors on `int()`, apply the NumPy 2 patch in README.md.
If anything fails on the day, show `demo_threshold72.png`, `dose_response_M1_vs_M3.png`
and `demo_fit.png` instead.

---

## Part 1: the structural claim (3 min, no fitting)

1. **What the two models are (20 s).**
   "Both are in PySB, same 8 species and 12 reactions. M1 is one reversible receptor.
   M3 is the same model plus a beta-arrestin2 brake on the occupied receptor. The only
   difference is `kcat_arr`: zero in M1, nonzero in M3."

2. **Run it.** `python demo.py`
   Read the two "reads as" columns out loud: M1 = down, down. M3 = down, then UP.
   Then the plot opens; point at the green curve crossing back over control.

3. **Say the caveat first.** "These rate constants are hand-picked. The claim here is structural:
   M1 cannot produce a sign flip at any parameter values, M3 can. The numbers come in part 2."

4. **Turn the knob.** The brake threshold is now `Km_arr`, in occupied receptors (of 100),
   not in dose units. That changed after his review.
   ```
   python demo.py --threshold 60     # flip moves to a lower dose
   python demo.py --threshold 85     # flip moves right, or disappears off the dose range
   ```
   "The threshold sets where the flip lands. Too low and the 150 dose stops looking protective,
   which contradicts the data. So the data bound this parameter from both sides."

5. **If he wants code.** `gev_models.py` is about 100 lines; the rule list is what he will read.
   `CHANGES_after_review.md` maps each of his comments to what changed.

---

## Part 2: the calibration (5 min)

6. **Start the live run first, talk while it goes.** `python demo_fit.py --live 150`
   That runs 3 chains x 150 iterations of PyDREAM on M3 live, about 70 seconds, so he sees
   real sampling rather than a saved number. While it runs, say what is being fitted:
   four targets, VCAM1 at 150 and 500, and the flip-protein mean at 150 and 500; Gaussian
   likelihood on the log2 scale; uniform priors in log10.
   Say plainly that 150 iterations is a demo burst, not the fit; it then loads the real
   1500-iteration run.

7. **The comparison table.** M1 max log-likelihood -27.2, AIC 66.4. M3 -10.2, AIC 40.4.
   Delta AIC = 26 for M3, after paying for 4 extra parameters.
   The line to say: "M1's best fit is a flat P, -0.31 at both doses. It never flips, at any
   parameter values. M3 gets -1.71 and +1.43 against data of -1.79 and +1.62."

8. **What three doses can and cannot pin down.** This is the honest part and the part that
   makes it systems biology rather than curve drawing.
   * PINNED: `fold_p`, `Km_arr`, `n_hill`.
   * unidentified: `KI_p`, `KA_p`, `kcat_arr`.
   "Three doses and ten parameters. The data pin the brake threshold and the flip amplitude.
   They say nothing about the arrestin catalytic rate. More doses would be needed for those."

9. **Say the convergence problem before he finds it.** Gelman-Rubin max is 1.61 on the M3 run;
   we want below 1.1. "This is 1500 iterations. 5000-plus before I quote any interval.
   The AIC gap is 26, so the model ranking is not in question, but the intervals are provisional."

10. **The posterior predictive figure.** `demo_fit.png` opens at the end: 40 posterior draws
    per model, best fit thick, data in black. Point at the right panel, where the grey M1 band
    is flat through zero and the green M3 band crosses it.

---

## Known misfits to raise yourself

* **VCAM1 at 500.** Data say -0.48, M3 says -0.91. The fit splits the difference between the
  two doses rather than reproducing the partial recovery. Either VCAM1 has its own recovery
  mechanism the model lacks, or the 500 point is noisy (FDR 0.41 in proteomics; the western
  blot only says "significantly lower", not how much).
* **Hill coefficient.** The posterior wants n of 7 to 15. That is steeper than a plain
  arrestin-recruitment step should give. Possibly receptor internalization and recycling,
  which we have not modeled.
* **Sigma is a placeholder.** 0.35 log2 units, not limma's per-protein standard errors.
  Those need exporting before any interval means anything.

## The questions for Dr. Harris
1. The Hill coefficient of 7 to 15 is doing a lot of work. Is adding internalization and
   recycling the right way to get that steepness from mechanism instead of from an exponent?
2. VCAM1 recovers partway at 500 in the data and the model does not reproduce it. Add a second
   VCAM1 input, or treat the 500 point as too noisy to fit given FDR 0.41?
3. Should `Arr_tot` be depletable rather than held constant? The QSSA assumes arrestin is in
   excess; at the top dose a large fraction of receptors is braked.

## If asked "why only five flip proteins"
Tnnc1 is excluded: troponin is absent from smooth muscle, so its 4.54 log2FC at 500 most likely
reflects cardiac tissue in the aortic arch samples. Histology is pending. Saa1 and Saa2 are
inverted-U but liver-derived acute-phase proteins, not artery-wall readouts. The five kept
(Nceh1, Comp, Serinc3, Aldh1l2, Sh3bgr) are all FDR below 1e-3 at both doses.
