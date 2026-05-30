"""
interactions.py — data layer for nuclear shell-model interactions.

Loads two-body matrix element (TBME) files (.int) from the CoSMo database
(A. Volya, github.com/alvolya/cosmo) and exposes them with PHASE-CORRECT
access, so the antisymmetry relation

    <(b a) J T| V |c d> = -(-1)^(j_a + j_b - J - T) <(a b) J T| V |c d>

is handled automatically. This is essential: the raw files store each TBME
once, in one particular orbit ordering, so naive lookups silently return 0
for the swapped ordering (which corrupts configuration mixing).

Common orbit convention used throughout: 1 = d3/2, 2 = d5/2, 3 = s1/2.
Files that use a different ordering (e.g. ab-initio JISP/N3LO) are remapped
to this convention on load via the CATALOG remap tables.
"""

import urllib.request

COSMO_RAW = "https://raw.githubusercontent.com/alvolya/cosmo/main/interactions"

# orbit index -> (name, n, l, j)   [sd shell, common convention]
SD_ORBITS = {1: ("d3/2", 0, 2, 1.5),
             2: ("d5/2", 0, 2, 2.5),
             3: ("s1/2", 1, 0, 0.5)}

# name -> (relative path in cosmo/interactions, remap to common convention or None)
CATALOG = {
    "USDB":       ("sd/usdb.int",       None),
    "USDA":       ("sd/usda.int",       None),
    "KUO":        ("sd/kuosd.int",      None),   # old G-matrix
    "WILDENTHAL": ("sd/w.int",          None),   # original USD
    "N3LO":       ("sdeff/N3LOA18.int", {1: 3, 2: 1, 3: 2}),  # chiral ab-initio
    "JISP":       ("sdeff/JISPA18.int", {1: 3, 3: 1, 5: 2}),  # ab-initio
}


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode()


def parse_int(text, remap=None):
    """Parse .int text -> (spe_dict, list of (a,b,c,d,J,T,V)).

    spe_dict maps orbit index -> single-particle energy (may be empty/zero for
    ab-initio files). Lines with fewer than 7 fields are skipped (handles the
    occasional malformed trailing line in some files).
    """
    started = False
    spe = {}
    rows = []
    for line in text.splitlines():
        if line.strip().startswith("!") or not line.strip():
            continue
        x = line.split()
        if not started:                       # first data line = header (nTBME + SPEs ...)
            started = True
            try:
                spe = {1: float(x[1]), 2: float(x[2]), 3: float(x[3])}
            except (IndexError, ValueError):
                spe = {}
            continue
        if len(x) < 7:
            continue
        try:
            a, b, c, d = int(x[0]), int(x[1]), int(x[2]), int(x[3])
            J, T, V = int(x[4]), int(x[5]), float(x[6])
        except ValueError:
            continue
        if remap:
            try:
                a, b, c, d = remap[a], remap[b], remap[c], remap[d]
            except KeyError:
                continue
        rows.append((a, b, c, d, J, T, V))
    return spe, rows


class Interaction:
    """A shell-model interaction: single-particle energies + TBME store."""

    def __init__(self, name, spe, rows, orbits=SD_ORBITS):
        self.name = name
        self.spe = spe
        self.orbits = orbits
        self._store = {}
        for a, b, c, d, J, T, V in rows:
            self._store[(a, b, c, d, J, T)] = V

    def j(self, a):
        return self.orbits[a][3]

    def getV(self, a, b, c, d, J, T):
        """Phase-correct <(a b) J T| V |(c d) J T>, for any orbit ordering."""
        def forms(p, q):
            ph = -((-1) ** int(self.j(p) + self.j(q) - J - T))
            return [((p, q), 1.0), ((q, p), ph)]
        for (ab, pab) in forms(a, b):
            for (cd, pcd) in forms(c, d):
                k = (ab[0], ab[1], cd[0], cd[1], J, T)
                if k in self._store:
                    return self._store[k] * pab * pcd
                k2 = (cd[0], cd[1], ab[0], ab[1], J, T)      # V is bra<->ket symmetric
                if k2 in self._store:
                    return self._store[k2] * pab * pcd
        return 0.0

    def monopole(self):
        """Monopole centroids M(a,b;T) = sum_J (2J+1) V(ab,ab;JT) / sum_J (2J+1)."""
        num, den = {}, {}
        for (a, b, c, d, J, T), V in self._store.items():
            if (a, b) == (c, d):
                key = (min(a, b), max(a, b), T)
                num[key] = num.get(key, 0.0) + (2 * J + 1) * V
                den[key] = den.get(key, 0.0) + (2 * J + 1)
        return {k: num[k] / den[k] for k in num}

    def canonical_tbme(self):
        """dict canonical (a<=b, c<=d, bra<=ket, J, T) -> raw V, for cross-interaction RMS.

        Uses the RAW stored value relabelled to a canonical key (NO antisymmetry
        phase). All CoSMo .int files share the same storage convention, so two
        interactions are compared in their common native ordering. (Forcing an
        index-sorted ordering with phases here would introduce convention-
        dependent sign flips that differ between a native file like USDB and a
        remapped one like JISP, and is unnecessary: the phase only matters when
        retrieving an arbitrary ordering inside the spectrum builder, via getV.)
        """
        out = {}
        for (a, b, c, d, J, T), V in self._store.items():
            br, ke = tuple(sorted((a, b))), tuple(sorted((c, d)))
            if br > ke:
                br, ke = ke, br
            out[(br[0], br[1], ke[0], ke[1], J, T)] = V
        return out


def load(name):
    """Download and parse a catalogued interaction by name (e.g. 'USDB', 'JISP')."""
    path, remap = CATALOG[name]
    spe, rows = parse_int(fetch(f"{COSMO_RAW}/{path}"), remap)
    return Interaction(name, spe, rows)


if __name__ == "__main__":
    import numpy as np
    usdb, usda = load("USDB"), load("USDA")
    print("USDB SPE:", {usdb.orbits[i][0]: round(usdb.spe[i], 3) for i in (1, 2, 3)})
    a = usdb.canonical_tbme(); b = usda.canonical_tbme()
    keys = [k for k in a if k in b]
    floor = np.sqrt(np.mean([(a[k] - b[k]) ** 2 for k in keys]))
    print(f"natural floor RMS(USDA,USDB) = {floor:.4f} MeV (expect 0.267)")
    print("phase check getV(d3/2,d5/2,d3/2,d5/2;J=1,T=0) =",
          round(usdb.getV(1, 2, 1, 2, 1, 0), 4), "(expect -6.0099)")
