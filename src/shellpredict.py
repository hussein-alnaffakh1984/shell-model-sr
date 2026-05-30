"""
shellpredict.py — a predictive & discovery toolkit for sd-shell effective
interactions, built on the validated pipeline of this repository.

What a researcher gets from it
------------------------------
Give it ANY interaction (a catalogued name, or a new .int file / URL) and it:

  1. PREDICTS the A = 18 spectra (18O, 18F) and compares with experiment & USDB;
  2. AUTO-FITS the per-isospin monopole correction Delta_T using our method
     (mean centroid difference to USDB) — no manual tuning;
  3. PREDICTS the monopole-corrected spectra (applies the DISCOVERED transferable
     T=0 correction);
  4. reports diagnostics: RMSE to USDB, monopole centroids, and where the
     interaction sits against the symbolic-regression information ceiling.

This is the workflow behind our paper, turned into a tool: drop in a new
ab-initio interaction and immediately see how it behaves and how the discovered
correction improves (or, honestly, sometimes worsens) the spectroscopy.

Everything reuses the validated modules interactions.py and spectrum.py, so the
numbers match the paper exactly (e.g. USDB 18O 2+ = 1.998 MeV).

Honest scope: the live calculator here is the validated two-valence-nucleon
(A = 18) engine. Extending prediction to arbitrary mid-shell nuclei needs a full
configuration-interaction engine — the next module on the roadmap.
"""

import copy
import numpy as np
import interactions as I
import spectrum as S

# ML information ceiling (cross-validated symbolic regression) and natural floor
ML_CEILING = 1.38
NATURAL_FLOOR = 0.267

# experimental A = 18 levels (MeV) used as the validation reference
EXP = {
    "18O": [(0.0, 0, 1), (1.98, 2, 1), (3.55, 4, 1)],
    "18F": [(0.0, 1, 0), (0.937, 3, 0), (1.042, 0, 1), (1.121, 5, 0)],
}


def load_any(spec, remap=None):
    """Load a catalogued interaction by name (e.g. 'USDB','JISP') or parse a new
    standard-ordered sd .int file from a path or URL."""
    if spec in I.CATALOG:
        return I.load(spec)
    text = I.fetch(spec) if spec.startswith("http") else open(spec).read()
    spe, rows = I.parse_int(text, remap)
    return I.Interaction(spec, spe, rows)


def rmse_to_usdb(inter, ref):
    a, b = inter.canonical_tbme(), ref.canonical_tbme()
    keys = [k for k in a if k in b]
    return float(np.sqrt(np.mean([(a[k] - b[k]) ** 2 for k in keys]))) if keys else None


def fit_correction(inter, ref):
    """Auto-fit per-isospin monopole correction Delta_T = <M_USDB - M_inter>_T.
    This reproduces the corrections derived by hand in the paper."""
    mi, mr = inter.monopole(), ref.monopole()
    out = {}
    for T in (0, 1):
        d = [mr[k] - mi[k] for k in mi if k[2] == T and k in mr]
        out[T] = float(np.mean(d)) if d else 0.0
    return out


def apply_correction(inter, delta):
    """Return a copy with each diagonal TBME of isospin T shifted by delta[T]."""
    c = copy.deepcopy(inter)
    c._store = {k: (v + delta.get(k[5], 0.0) if (k[0], k[1]) == (k[2], k[3]) else v)
                for k, v in c._store.items()}
    return c


def predict_A18(inter, spe):
    """Return {'18O': [(E,J,T)...], '18F': [...]} for this interaction."""
    return {"18O": S.o18(inter, spe=spe)[1], "18F": S.f18(inter, spe=spe)[1]}


def _match(levels, J, T):
    for e, j, t in levels:
        if j == J and t == T:
            return e
    return None


def report(spec, use_usdb_spe=None):
    """Print a full prediction + discovery report for an interaction.

    `spec` is a catalogue name, path, or URL. `use_usdb_spe` overrides the
    single-particle energies with USDB's (default: True for ab-initio files
    whose stored SPEs are zero, False for phenomenological fits)."""
    usdb = I.load("USDB")
    inter = load_any(spec)
    name = inter.name
    ab_initio = all(abs(inter.spe.get(i, 0.0)) < 1e-9 for i in (1, 2, 3))
    if use_usdb_spe is None:
        use_usdb_spe = ab_initio
    spe = usdb.spe if use_usdb_spe else inter.spe

    rmse = rmse_to_usdb(inter, usdb)
    delta = fit_correction(inter, usdb)
    base = predict_A18(inter, spe)
    corr = predict_A18(apply_correction(inter, delta), spe)

    print(f"\n{'='*60}\n  PREDICTION REPORT — {name}\n{'='*60}")
    if rmse is not None:
        pos = "below floor" if rmse <= NATURAL_FLOOR else (
              "near floor" if rmse < 0.5 else
              "mid-range" if rmse < ML_CEILING else "above ML ceiling")
        print(f"  RMSE to USDB           : {rmse:.3f} MeV   [{pos}]")
        print(f"  (natural floor {NATURAL_FLOOR}, ML ceiling ~{ML_CEILING})")
    print(f"  auto-fitted monopole Δ : T=0 {delta[0]:+.3f},  T=1 {delta[1]:+.3f} MeV")
    print(f"  single-particle energies: {'USDB (ab-initio file)' if use_usdb_spe else name}")

    for nuc in ("18O", "18F"):
        print(f"\n  {nuc}:  Jπ(T)      exp     predict   +corr")
        seen = set()
        for e_exp, J, T in EXP[nuc]:
            b = _match(base[nuc], J, T)
            c = _match(corr[nuc], J, T)
            tag = (J, T)
            if tag in seen:
                continue
            seen.add(tag)
            bs = f"{b:6.2f}" if b is not None else "   -- "
            cs = f"{c:6.2f}" if c is not None else "   -- "
            print(f"     {J}+(T{T})   {e_exp:6.2f}   {bs}   {cs}")
    print()
    return {"name": name, "rmse": rmse, "delta": delta, "base": base, "corrected": corr}


if __name__ == "__main__":
    import sys
    targets = sys.argv[1:] or ["USDB", "JISP", "N3LO"]
    for tgt in targets:
        report(tgt)
