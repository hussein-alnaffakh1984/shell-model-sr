"""
features.py — feature engineering for symbolic regression on TBME.

Turns an Interaction into a tidy table: one row per (canonical) two-body
matrix element, with its value V and a set of input features. Two feature
sets are exposed:

  QUANTUM_FEATURES : raw quantum numbers only (J, T, and the (n, l, j) of the
                     four single-particle orbits). Symbolic regression on
                     these alone plateaus at the "information ceiling".

  PHYSICS_FEATURES : the above plus physically motivated quantities --
                     pairing flag, diagonal flag, monopole centroid of the
                     (a,b) pair, spin-orbit label (j - l), coupling features.
                     The monopole centroid in particular carries most of the
                     learnable signal.
"""

import numpy as np
import pandas as pd

QUANTUM_FEATURES = ["J", "T", "ja", "jb", "jc", "jd",
                    "la", "lb", "lc", "ld", "na", "nb", "nc", "nd"]

PHYSICS_FEATURES = QUANTUM_FEATURES + [
    "twoJ1", "JJ1", "pairing", "diagonal", "same_bra", "same_ket",
    "stretch", "sumj", "dj_bra", "dj_ket",
    "ls_a", "ls_b", "ls_c", "ls_d", "centroid_ab", "spe_avg",
]


def build_features(inter):
    """Return a pandas DataFrame of TBME + features for one interaction."""
    cen = inter.monopole()
    spe = inter.spe
    orb = inter.orbits
    rows = []
    for (a, b, c, d, J, T), V in inter.canonical_tbme().items():
        (_, na, la, ja) = orb[a]
        (_, nb, lb, jb) = orb[b]
        (_, nc, lc, jc) = orb[c]
        (_, nd, ld, jd) = orb[d]
        diagonal = int((a, b) == (c, d))
        rows.append(dict(
            interaction=inter.name, a=a, b=b, c=c, d=d, J=J, T=T, V=V,
            ja=ja, jb=jb, jc=jc, jd=jd, la=la, lb=lb, lc=lc, ld=ld,
            na=na, nb=nb, nc=nc, nd=nd,
            twoJ1=2 * J + 1, JJ1=J * (J + 1),
            pairing=int(J == 0 and T == 1 and diagonal),
            diagonal=diagonal, same_bra=int(a == b), same_ket=int(c == d),
            stretch=int(abs((ja + jb) - J) < 1e-6),
            sumj=ja + jb + jc + jd, dj_bra=abs(ja - jb), dj_ket=abs(jc - jd),
            ls_a=ja - la, ls_b=jb - lb, ls_c=jc - lc, ls_d=jd - ld,
            centroid_ab=cen.get((min(a, b), max(a, b), T), 0.0),
            spe_avg=0.25 * (spe[a] + spe[b] + spe[c] + spe[d]) if spe else 0.0,
        ))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    from interactions import load
    df = build_features(load("USDB"))
    print("USDB feature table:", df.shape, "(expect 63 rows)")
    print("feature columns:", len(PHYSICS_FEATURES))
    print(df[["a", "b", "c", "d", "J", "T", "V", "centroid_ab", "pairing"]].head(6).to_string(index=False))
