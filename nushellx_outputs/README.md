# NuShellX@MSU output files

This folder holds the spectra behind Table 2 and Supplementary Table S1.

## Status

| Nucleus | File(s) | Source |
|---------|---------|--------|
| 17O | `17O_USDA.lpt`, `17O_N3LOstar.lpt` | **Real raw states**, formatted from the team's Extended Generalization Report (Tables 10–11). |
| 19F | `19F_USDA.lpt`, `19F_N3LOstar.lpt` | **Real raw states**, formatted from the team's Extended Generalization Report (Tables 12–13). |
| 18O | `18O_levels.txt` | Level energies only (from comparison tables). Replace with full `.lpt`. |
| 20Ne | `20Ne_levels.txt` | Level energies only. Replace with full `.lpt`. |
| 24Mg | `24Mg_levels.txt` | Level energies only. Replace with full `.lpt`. |

## What the team still needs to provide

To make the reproduction package fully self-contained, replace the three
`*_levels.txt` files with the complete NuShellX `.lpt` output, named:

```
18O_USDA.lpt   18O_N3LOstar.lpt
20Ne_USDA.lpt  20Ne_N3LOstar.lpt
24Mg_USDA.lpt  24Mg_N3LOstar.lpt
```

`N3LOstar` = chiral N3LO with USDA single-particle energies and the
18O-tuned correction (written "N3LO*" in the paper).

## Honesty note

The 17O and 19F `.lpt` files are **reformatted** from the team's tabulated
output (the raw numbers are exactly as reported). They are faithful to the
reported values but are not the original binary NuShellX files; the team
retains those and they are the authoritative source for any referee check.
