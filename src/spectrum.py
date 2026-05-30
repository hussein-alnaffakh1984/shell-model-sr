"""
spectrum.py — two-particle shell-model spectra for A = 18 nuclei.

For two valence nucleons (A = 18 on a 16O core) the configuration-interaction
problem is small enough to diagonalise exactly with NumPy. For each total
angular momentum J and isospin T we build the Hamiltonian in the basis of
normalised antisymmetrised pairs |(a<=b) J T>,

    H[(a,b),(c,d)] = <(a b) J T| V |(c d) J T> + (eps_a + eps_b) delta,

using the PHASE-CORRECT matrix elements from interactions.Interaction.getV,
then diagonalise. This is validated against experiment:

    18O (T=1): 2+ = 1.998 MeV (exp 1.98), 4+ = 3.527 MeV (exp 3.55)
    18F (T=0,1): 1+ ground state (correct), 3+ = 0.945 MeV (exp 0.937)

Selection rules:
    same orbit (a==b): allowed only if (J + T) is odd   [antisymmetry]
    different orbits : |j_a - j_b| <= J <= j_a + j_b
"""

import numpy as np


def _pairs(orbs):
    return [(a, b) for i, a in enumerate(orbs) for b in orbs[i:]]


def channel(inter, J, T, spe=None):
    """Eigenvalues (sorted) of the (J, T) two-particle Hamiltonian.

    `spe` optionally overrides the interaction's single-particle energies
    (useful for ab-initio files whose stored SPEs are zero)."""
    spe = inter.spe if spe is None else spe
    orbs = sorted(inter.orbits)
    basis = []
    for (a, b) in _pairs(orbs):
        ja, jb = inter.j(a), inter.j(b)
        if a == b:
            if (J + T) % 2 == 0 or J > 2 * ja:
                continue
        else:
            if J < abs(ja - jb) or J > ja + jb:
                continue
        basis.append((a, b))
    if not basis:
        return None
    n = len(basis)
    H = np.zeros((n, n))
    for i, (a, b) in enumerate(basis):
        for k, (c, d) in enumerate(basis):
            H[i, k] = inter.getV(a, b, c, d, J, T)
        H[i, i] += spe[a] + spe[b]
    return np.sort(np.linalg.eigvalsh(H))


def spectrum(inter, channels, spe=None, jmax=6):
    """Return {(J, T): sorted eigenvalues} over the requested isospin channels."""
    out = {}
    for T in channels:
        for J in range(0, jmax):
            ev = channel(inter, J, T, spe=spe)
            if ev is not None:
                out[(J, T)] = ev
    return out


def levels(inter, channels, spe=None, emax=6.0):
    """Excitation energies as a list of (E_star, J, T), ground state at 0."""
    sp = spectrum(inter, channels, spe=spe)
    gs = min(ev[0] for ev in sp.values())
    rows = []
    for (J, T), ev in sp.items():
        for e in ev:
            ex = e - gs
            if -1e-6 <= ex <= emax:
                rows.append((round(ex, 3), J, T))
    return gs, sorted(rows)


def o18(inter, spe=None):
    """18O: two neutrons (T = 1 only)."""
    return levels(inter, [1], spe=spe)


def f18(inter, spe=None):
    """18F: proton + neutron (T = 0 and T = 1)."""
    return levels(inter, [0, 1], spe=spe)


if __name__ == "__main__":
    from interactions import load
    usdb = load("USDB")

    gs, lv = o18(usdb)
    print("USDB 18O  (exp: 2+ 1.98, 4+ 3.55)")
    for ex, J, T in lv:
        if ex < 4:
            print(f"   {J}+  E* = {ex:.3f} MeV")

    gs, lv = f18(usdb)
    print("USDB 18F  (exp: 1+ gs, 3+ 0.94, 0+(T1) 1.04)")
    for ex, J, T in lv:
        if ex < 1.6:
            print(f"   J={J} T={T}  E* = {ex:.3f} MeV")
