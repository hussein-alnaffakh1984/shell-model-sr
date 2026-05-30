# Shell-Model Effective Interaction: ML Ceiling, *ab-initio* Baseline, and Spectra

A reproducible study of **how much of the sd-shell effective interaction is
"learnable"** from quantum-number features, and where the missing physics comes
from. We combine **symbolic regression** (the ML core) with **modern ab-initio
interactions** and validate everything against real nuclear spectra.

All data come from the [CoSMo](https://github.com/alvolya/cosmo) interaction
database (A. Volya). Every number below is reproduced by `notebooks/run_all.ipynb`.

## Key results

| quantity | value | meaning |
|---|---|---|
| natural floor RMS(USDA, USDB) | **0.267 MeV** | intrinsic calibration scatter |
| ML ceiling (symbolic regression, CV) | **~1.4 MeV** | best reachable with quantum-number features |
| discovered leading form | **V ≈ diagonal × monopole-centroid** | what SR converges to |
| best ab-initio baseline (JISP, 0 params) | **0.466 MeV** | far below the ML ceiling |
| old G-matrix (Kuo) baseline | 1.081 MeV | a poor baseline |
| JISP + 2-parameter monopole correction | **0.404 MeV** | T=0 part transfers sd→fp |
| USDB ¹⁸O 2⁺ / ¹⁸F 1⁺ g.s. | 1.998 MeV / correct | spectrum machinery validated |

**Headline:** symbolic regression *measures the information ceiling* of
quantum-number features; modern ab-initio supplies the missing physics
(reproducing sd spectra with **zero** fitting); a naive monopole RMS-correction
reduces TBME error but does **not** improve—and can worsen—spectroscopic
observables. This last point is an honest, cautionary result about the
TBME-RMS metric.

## Repository layout

```
shell-model-sr/
├── README.md
├── PAPER_OUTLINE.md         # planned paper: story, sections, figures
├── src/
│   ├── interactions.py      # data layer: load/parse .int, phase-correct TBME, monopole
│   ├── features.py          # feature engineering for symbolic regression
│   ├── symbolic.py          # symbolic-regression pipeline (PySR) — the ML core
│   └── spectrum.py          # validated A=18 (¹⁸O, ¹⁸F) shell-model spectra
└── notebooks/
    └── run_all.ipynb        # runs the whole pipeline, prints every verified number
```

## Quick start

```bash
pip install numpy pandas scikit-learn pysr   # PySR installs Julia on first import
```

```python
import sys; sys.path.insert(0, "src")
import interactions as I
from features import build_features, PHYSICS_FEATURES
from symbolic import fit_sr, cv_rmse
import spectrum as S

usdb = I.load("USDB")
print(S.o18(usdb))                       # 18O spectrum
df = build_features(usdb)
print(cv_rmse(df, PHYSICS_FEATURES))     # ML ceiling ~1.4 MeV
```

Or just open **`notebooks/run_all.ipynb`** on Kaggle (Internet = ON) and *Run All*.

## Notes

* `interactions.canonical_tbme()` compares interactions in their common native
  ordering (raw values); `interactions.getV()` applies the antisymmetry phase
  and is used only inside the spectrum builder. Mixing the two conventions is a
  classic source of silent sign errors — see the docstrings.
* PySR is configured for speed and clean output: `parallelism="multithreading"`,
  `verbosity=0`, `progress=False`.

## Data & references

* Data: CoSMo, A. Volya — https://github.com/alvolya/cosmo
* ab-initio sd interactions: Dikmen *et al.*, Phys. Rev. C **91**, 064301 (2015)
* Related (UQ refit of USDB): Gorton & Kravvaris, arXiv:2503.11889 (2025)
