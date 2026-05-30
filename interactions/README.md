# Interaction files (.int)

Standard sd-shell format. Orbit order: 1 = 0d3/2, 2 = 0d5/2, 3 = 1s1/2.
First data line = header: (-)count, three single-particle energies, A_core, A_ref, mass-exponent.
Mass dependence V(A) = V(18) * (18/(16+n))^0.30 applied by NuShellX.

- USDA.int            phenomenological reference (Brown & Richter 2006)
- N3LO_usda_spe.int   N3LO chiral NN (Entem-Machleidt) two-body + USDA single-particle energies
- N3LO_plus_corr.int  N3LO + transferable per-isospin monopole correction
- N3LO_optimized.int  N3LO after 6-channel tuning on 18O (generalization study)
