# Paper outline (working draft)

> Honest plan for the paper, based on what the analysis actually showed.
> Keep claims matched to the verified numbers — no overselling.

## Working title

*"How much of the shell-model effective interaction is learnable? An
information-ceiling analysis with symbolic regression and ab-initio baselines."*

(Alternative: *"Symbolic regression meets ab initio: the information content of
the sd-shell effective interaction."*)

## One-sentence thesis

Symbolic regression quantifies an **information ceiling** for describing the
sd-shell effective interaction from quantum-number features; the physics beyond
that ceiling is supplied, parameter-free, by modern ab-initio interactions — and
reducing two-body-matrix-element (TBME) RMS by monopole corrections does not, by
itself, improve spectroscopic observables.

## Abstract (sketch)

We ask how much of the empirical sd-shell interaction (USDB) is determined by the
quantum numbers of the interacting orbits. Using symbolic regression we find an
**information ceiling** of ~1.4 MeV (cross-validated), well above the 0.267 MeV
natural calibration floor, with the leading learnable structure being
`V ≈ diagonal × monopole-centroid`. We then show that modern chiral ab-initio
interactions (JISP, N3LO) reproduce USDB two-body matrix elements to ~0.47 MeV
and the ¹⁸O / ¹⁸F spectra **without any fitting**, whereas an older G-matrix
(Kuo) does not. The empirical residual is monopole-dominated; its **isoscalar
(T=0)** component (~−0.4 MeV) is approximately transferable from the sd to the fp
shell, while the isovector component is not. Finally, we demonstrate that a
monopole correction that minimizes TBME RMS does **not** improve — and can
degrade — spectroscopic observables (e.g. the ¹⁸F isobaric-analog state),
a cautionary result for data-driven interaction fitting.

## Three results (the spine)

1. **The ML ceiling.** SR on quantum-number features plateaus at CV-RMSE ~1.4 MeV
   (in-sample ~1.1); discovered form `diagonal × centroid_ab`. Natural floor
   (USDA↔USDB) = 0.267 MeV. → quantum numbers alone cannot reach the floor.

2. **ab-initio is the missing physics.** Zero-parameter baseline scan vs USDB:
   Kuo 1.081, Wildenthal 1.020, N3LO 0.484, **JISP 0.466**. Bare JISP reproduces
   ¹⁸O (2⁺ ≈ 1.7–2.0) and ¹⁸F (1⁺ g.s.; 3⁺ ≈ 0.99) spectra unfitted.

3. **Monopole correction: transferable but not spectroscopically helpful.**
   JISP + 2 per-T constants → 0.404 MeV. Shifts: T0 = −0.42, T1 = +0.34 (sd).
   Transfer test sd→fp: T0 transfers (−0.42 → −0.33), T1 does not (+0.34 → −0.06;
   caveat: fp baseline kb3g is circular with GXPF1A). The correction leaves
   ¹⁸F internal excitations unchanged and **worsens** the IAS (1.6 → 2.4 vs exp
   1.04) — minimizing TBME RMS ≠ improving observables.

## Section plan

1. **Introduction** — effective interactions; "fit the numbers" (e.g. USD,
   USDB, and the UQ refit of Gorton & Kravvaris 2025) vs "find the structure";
   why an information-content question is worth asking.
2. **Data and methods** — CoSMo database; feature engineering; symbolic
   regression (PySR) setup and cross-validation; zero-parameter baseline
   comparison; monopole decomposition; A=18 configuration-interaction spectra
   (phase-correct TBME).
3. **Results**
   - 3.1 Information ceiling and discovered form.
   - 3.2 ab-initio baselines and parameter-free spectra.
   - 3.3 Monopole correction and sd→fp transferability.
   - 3.4 Spectroscopic test: TBME-RMS vs observables.
4. **Discussion** — what the ceiling means; why the TBME metric misleads;
   implications for data-driven / ML-assisted interaction building; outlook for
   predictions toward the dripline.
5. **Conclusion.**

## Figures (planned)

- **F1** SR Pareto front + ceiling (1.4) vs floor (0.267) bands.
- **F2** baseline scan bar chart (Kuo, Wildenthal, N3LO, JISP) with the ceiling/floor lines.
- **F3** monopole correction per (orbit pair, T), sd vs fp (shows T=0 transfer).
- **F4** ¹⁸O and ¹⁸F level schemes: experiment vs USDB vs JISP vs JISP+correction.

## Honest caveats to state explicitly

- fp transferability is indicative, not definitive (kb3g↔GXPF1A circularity; no
  independent chiral fp interaction in CoSMo).
- 2-parameter monopole correction reaches 0.40 MeV, still above the 0.267 floor —
  an irreducible residual remains.
- The monopole correction's failure to improve spectra is a real (publishable)
  limitation, not a tuning issue.

## Open / future

- Independent ab-initio fp interaction → clean isovector transfer test.
- Correction fit to **observables** rather than TBME RMS.
- Uncertainty propagation (USDBUQ-style σ as SR weights).
