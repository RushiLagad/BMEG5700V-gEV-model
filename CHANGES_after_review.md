# What changed after Dr. Harris's review (Nov 5)

## His comment 1: the binding rule was a synthesis rule
v1:  `gEV + CaSR >> gEV + CaSR%gEV`
     gEV appears on both sides, so this *creates* a bound complex and never consumes ligand.
     His note: "synthesis + binding". Not a binding rule.
v2:  `gEV + CaSR | gEV%CaSR`  with kf, kr   (`|` = bidirectional in PySB; `>>` = unidirectional)
     The dose is clamped the proper way: `Initial(gEV(b=None), gEV_0, fixed=True)`.

## His comment 2: the brake rate mixed a Hill term with mass action
v1 rate on `occ + Arr >> occ%Arr`:
     k1 * [gEV]^n / (K^n + [gEV]^n)   times   [occ] * [Arr]   (PySB adds the reactants)
     A Hill function in a species that is not a reactant, times mass action in two that are.
     He expanded it on page 2 to show the full third-order rate, then derived Michaelis-Menten
     on pages 2-3 and boxed the Hill form:  S -> P with keff = kcat * Etot * [S]^n / (Km^n + [S]^n).
v2:  arrestin is the enzyme in excess (Arr_tot, no longer an explicit species);
     the OCCUPIED RECEPTOR is the substrate; the rule is unimolecular:
         CaSR(b=1,arr=off)%gEV >> CaSR(b=1,arr=on)%gEV
         flux = kcat_arr * Arr_tot * occ^n / (Km_arr^n + occ^n)
     PySB multiplies a unimolecular rule's rate by [occ], so the Expression is flux / occ.
     The reset rule now releases the ligand too (his rule 3): braked -> CaSR + gEV.

## Consequence
The "positive feedback" I worried about in v1 was an artifact of keying the switch to a
species that braking does not deplete. Keyed to the occupied receptor, braking consumes its
own substrate, so it is self-limiting. The flip survives, but in a narrower parameter region:
     Km_arr ~ 72 (occupied receptors), kcat_arr 0.3, n 12, kreset 0.5, KA_p 10
     -> P: 0.57x at 150, 1.81x at 500.   VCAM1 low at both doses.   M1 still cannot flip.

## Still to confirm with him
* n = 12 is steep. Is a Hill coefficient that high defensible for arrestin recruitment, or is it
  a sign that another mechanism (e.g. receptor internalisation with slow recycling) should carry
  part of the switch?
* Whether Arr_tot should stay a constant or come back as a depletable species.
