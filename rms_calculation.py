#!/usr/bin/env python3
"""
rms_calculation.py  --  regenerate Table 2 of the manuscript from NuShellX output.

RMS definition (matches the Methods section):

    RMS = sqrt( (1/N) * sum_i (Ex_i^calc - Ex_i^exp)^2 )

  * Ex are EXCITATION energies (relative to the computed ground state), not
    absolute binding energies.
  * The sum runs over the explicitly listed yrast states in SELECTED_STATES.
  * The ground state contributes a zero residual (INCLUDE_GROUND_STATE).
  * Non-sd intruder states (e.g. 17O 3/2+) are excluded (EXCLUDE_INTRUDERS).

For 17O and 19F the repository ships real .lpt files (formatted from the team's
NuShellX output); this script reads them directly and reproduces the published
rms. For 18O/20Ne/24Mg the repository currently ships *_levels.txt summaries;
the script falls back to those until full .lpt files are added.

Usage:
    python rms_calculation.py
"""

import math, os, re, glob

INCLUDE_GROUND_STATE = True
EXCLUDE_INTRUDERS    = True
DIR = os.path.join(os.path.dirname(__file__), "nushellx_outputs")

# experimental excitation energies (MeV); sources cited in the paper
EXP = {
 "17O": {"5/2+":0.000,"1/2+":0.871,"3/2+":5.085},
 "18O": {"0+":0.000,"2+":1.982,"4+":3.555,"0+_2":3.634,"2+_2":3.920,"3+":5.378},
 "19F": {"1/2+":0.000,"5/2+":0.197,"3/2+":1.346,"9/2+":2.780,"7/2+":4.033},
 "20Ne":{"0+":0.000,"2+":1.634,"4+":4.248},
 "24Mg":{"0+":0.000,"2+":1.369,"4+":4.123},
}
# yrast states entering the rms (the single source of truth)
SELECTED = {
 "17O": ["1/2+"],                       # 3/2+ excluded (p-sd intruder)
 "18O": ["2+","4+"],
 "19F": ["5/2+","3/2+","9/2+","7/2+"],
 "20Ne":["2+","4+"],
 "24Mg":["2+","4+"],
}
INTRUDERS = {("17O","3/2+")}

def jlabel(twoJ, parity):
    p = "+" if int(parity) > 0 else "-"
    twoJ = int(twoJ)
    return (f"{twoJ//2}" if twoJ % 2 == 0 else f"{twoJ}/2") + p

def parse_lpt(path):
    """Read a reformatted NuShellX .lpt: columns N NJ E Ex 2J 2T parity file."""
    out = []
    for line in open(path):
        if line.strip().startswith("!") or not line.strip():
            continue
        s = line.split()
        if len(s) < 7:
            continue
        try:
            ex = float(s[3]); twoJ = int(s[4]); par = s[6]
        except ValueError:
            continue
        par_val = +1 if par in ("+","1","+1") else -1
        out.append({"Ex": ex, "J": jlabel(twoJ, par_val)})
    return out

def lowest_ex(states, label):
    cands = [st["Ex"] for st in states if st["J"] == label]
    return min(cands) if cands else None

def rms(res):
    return math.sqrt(sum(r*r for r in res)/len(res)) if res else float("nan")

def calc_from_lpt(nuc, tag):
    path = os.path.join(DIR, f"{nuc}_{tag}.lpt")
    if not os.path.exists(path):
        return None
    states = parse_lpt(path)
    res = [0.0] if INCLUDE_GROUND_STATE else []
    for lab in SELECTED[nuc]:
        if EXCLUDE_INTRUDERS and (nuc, lab) in INTRUDERS:
            continue
        calc = lowest_ex(states, lab)
        if calc is None:
            return None
        res.append(calc - EXP[nuc][lab])
    return rms(res)

def calc_from_levels(nuc, col):
    """Fallback: read *_levels.txt (Jpi Exp USDA N3LO*)."""
    path = os.path.join(DIR, f"{nuc}_levels.txt")
    if not os.path.exists(path):
        return None
    vals = {}
    for line in open(path):
        if line.strip().startswith("#") or not line.strip():
            continue
        s = line.split()
        lab = s[0].replace("_1","")
        vals[lab] = {"USDA": float(s[2]), "N3LO*": float(s[3])}
    res = [0.0] if INCLUDE_GROUND_STATE else []
    for lab in SELECTED[nuc]:
        key = lab if lab in vals else lab+"_1"
        if key not in vals:
            return None
        res.append(vals[key][col] - EXP[nuc][lab])
    return rms(res)

def main():
    print("RMS = sqrt( (1/N) Σ (Ex_calc - Ex_exp)^2 ),  excitation energies only,")
    print(f"INCLUDE_GROUND_STATE={INCLUDE_GROUND_STATE}, EXCLUDE_INTRUDERS={EXCLUDE_INTRUDERS}\n")
    print(f"{'Nucleus':8}{'USDA':>10}{'N3LO*':>10}   source")
    for nuc in ["17O","18O","19F","20Ne","24Mg"]:
        ru = calc_from_lpt(nuc,"USDA");  rn = calc_from_lpt(nuc,"N3LOstar")
        src = ".lpt"
        if ru is None or rn is None:
            ru = calc_from_levels(nuc,"USDA"); rn = calc_from_levels(nuc,"N3LO*")
            src = "levels.txt (awaiting full .lpt)"
        su = f"{ru:.3f}" if ru is not None else "  n/a"
        sn = f"{rn:.3f}" if rn is not None else "  n/a"
        print(f"{nuc:8}{su:>10}{sn:>10}   {src}")

if __name__ == "__main__":
    main()
