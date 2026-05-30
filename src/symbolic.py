"""
symbolic.py — symbolic-regression pipeline (the ML core).

Wraps PySR with the settings that were validated on Kaggle:

  * parallelism="multithreading"  -> uses all CPU cores (fast)
  * verbosity=0, progress=False   -> suppresses the Julia progress table that
                                     triggers the cosmetic UnicodeDecodeError
  * optional per-element weights  -> hook for USDBUQ-style uncertainties
                                     (weight = 1/sigma^2)

Typical results on USDB (sd shell):
  in-sample RMSE   ~ 1.3 MeV
  CV-RMSE          ~ 1.4 MeV   <- the "information ceiling" of quantum-number
                                  features; the discovered form is essentially
                                  V ~ diagonal * monopole_centroid.

PySR needs a working Julia backend (installed automatically on first import,
internet required). This module imports PySR lazily so the rest of the repo
can be used without it.
"""

import numpy as np

BINARY_OPS = ["+", "-", "*", "/"]
UNARY_OPS = ["square", "cube", "abs"]
WEIGHTED_LOSS = "loss(pred, y, w) = w * (pred - y)^2"


def _make_model(niterations, populations, maxsize, weighted):
    from pysr import PySRRegressor
    kw = dict(
        niterations=niterations, populations=populations, population_size=35,
        maxsize=maxsize, binary_operators=BINARY_OPS, unary_operators=UNARY_OPS,
        model_selection="best", progress=False, verbosity=0,
        parallelism="multithreading",
    )
    if weighted:
        kw["elementwise_loss"] = WEIGHTED_LOSS
    return PySRRegressor(**kw)


def _fit(model, X, y, w, features):
    if w is not None:
        model.fit(X, y, weights=w, variable_names=list(features))
    else:
        model.fit(X, y, variable_names=list(features))
    return model


def fit_sr(df, features, target="V", weights=None,
           niterations=60, populations=20, maxsize=26):
    """Fit a symbolic-regression model; return (model, in_sample_rmse)."""
    X = df[features].values.astype(float)
    y = df[target].values.astype(float)
    w = None if weights is None else np.asarray(weights, float)
    model = _fit(_make_model(niterations, populations, maxsize, w is not None),
                 X, y, w, features)
    rmse = float(np.sqrt(np.mean((model.predict(X) - y) ** 2)))
    return model, rmse


def cv_rmse(df, features, target="V", weights=None, k=3,
            niterations=35, populations=15, maxsize=22, seed=1):
    """Cross-validated RMSE — the honest generalization estimate."""
    from sklearn.model_selection import KFold
    X = df[features].values.astype(float)
    y = df[target].values.astype(float)
    w = None if weights is None else np.asarray(weights, float)
    errs = []
    for tr, te in KFold(k, shuffle=True, random_state=seed).split(X):
        m = _make_model(niterations, populations, maxsize, w is not None)
        wt = None if w is None else w[tr]
        m = _fit(m, X[tr], y[tr], wt, features)
        errs.append(np.sqrt(np.mean((m.predict(X[te]) - y[te]) ** 2)))
    return float(np.mean(errs)), float(np.std(errs))


def best_equations(model, n=6):
    """Return the Pareto front (complexity / loss / equation) as a DataFrame."""
    return model.equations_[["complexity", "loss", "equation"]].head(n)


if __name__ == "__main__":
    from interactions import load
    from features import build_features, PHYSICS_FEATURES
    df = build_features(load("USDB"))
    print("Fitting SR on USDB V (physics features)...")
    model, rmse = fit_sr(df, PHYSICS_FEATURES)
    print(f"in-sample RMSE = {rmse:.3f} MeV")
    print(best_equations(model).to_string(index=False))
    mu, sd = cv_rmse(df, PHYSICS_FEATURES)
    print(f"CV-RMSE = {mu:.3f} +/- {sd:.3f} MeV  (expect ~1.4 = information ceiling)")
