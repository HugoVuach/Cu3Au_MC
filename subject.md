# Chapitre 1

# Study of the order-disorder transition of the Cu₃Au compound by the Monte Carlo method on a rigid lattice (FD) \*\*

## Objective

You will here **write yourselves** a Metropolis Monte Carlo code on a rigid lattice, based on a **simplified pair model** (a two-species generalized Ising-type model), restricted to first and second neighbors.

## 1.1 Physical model

The Cu₃Au compound crystallizes, in its ordered state (L1₂ structure), in a cubic cell of parameter \(a\) containing 4 sites:

- **Au**: vertex \((0,0,0)\)
- **Cu**: face centers \(\left(0,\frac12,\frac12\right)\), \(\left(\frac12,0,\frac12\right)\), \(\left(\frac12,\frac12,0\right)\)

The underlying lattice (site positions, independent of the species occupying them) is a **face-centered cubic (fcc)** lattice. The order-disorder transition consists in randomly permuting Au and Cu atoms on this *fixed* fcc lattice (frozen positions – hence “rigid lattice”); only the occupation (Au or Cu) of each site varies.

Build a cubic box of \(6 \times 6 \times 6\) cells (in units of \(a\)), i.e. \(6^3 \times 4 = 864\) sites in total, with periodic boundary conditions (PBC) in the three directions.

At \(T = 0\ \mathrm{K}\), the reference configuration is the ordered L1₂ structure described above (one quarter of the sites are Au, three quarters are Cu, distributed exactly as in the elementary cell repeated \(6 \times 6 \times 6\) times).

For a given fcc site, we consider two neighbor shells:

| Shell | Distance | Number of neighbors \(z\) |
|---|---|---|
| 1st neighbors | \(a/\sqrt{2}\) | 12 |
| 2nd neighbors | \(a\) | 6 |

*Practical advice*: the simplest approach is to generate, once and for all at the beginning of the program, the list of the 12 indices of the 1st neighbors and of the 6 indices of the 2nd neighbors for each site, applying PBC (modulo 6 on each cell coordinate). This neighbor table never changes during the simulation – only the occupation of the sites changes.

The total energy of the system is a sum over pairs of first- and second-neighbor sites, each pair contributing according to the nature of the two atoms occupying it:

$$
E =
\sum_{\langle i,j\rangle_1} J_{\sigma_i\sigma_j}
+
\sum_{\langle i,j\rangle_2} I_{\sigma_i\sigma_j}
\tag{1.1}
$$

where \(\sigma_i \in \{\mathrm{Au},\mathrm{Cu}\}\) is the species occupying site \(i\), the first sum runs over pairs of 1st neighbors, the second over pairs of 2nd neighbors, and

$$
J_{\sigma_i\sigma_j}
\in
\{J_{\mathrm{AuAu}},J_{\mathrm{CuCu}},J_{\mathrm{AuCu}}\},
\qquad
I_{\sigma_i\sigma_j}
\in
\{I_{\mathrm{AuAu}},I_{\mathrm{CuCu}},I_{\mathrm{AuCu}}\}
\tag{1.2}
$$

according to the nature of the two atoms of the pair.

| Parameter | Value (eV) |
|---|---:|
| \(J_{\mathrm{AuAu}}\) | −0.635 |
| \(J_{\mathrm{CuCu}}\) | −0.582 |
| \(J_{\mathrm{AuCu}}\) | −0.689 |
| \(I_{\mathrm{AuAu}}\) | 0.000 |
| \(I_{\mathrm{CuCu}}\) | 0.000 |
| \(I_{\mathrm{AuCu}}\) | +0.0063 |

## 1.2 Simulation of disorder

1. Write the Monte Carlo program on a rigid lattice adapted to the study of the order-disorder transition of Cu₃Au: construction of the crystal lattice and of the neighbor table (1st and 2nd neighbors with PBC), potential energy of the system from the pair model above, exchange of atoms, Metropolis algorithm, …

2. Write the subprogram computing the Cowley short-range order parameters at first and second neighbors (function `cowley`).

   Reminder: the Cowley parameter for shell \(n\) is written, for a reference site of species Au, as

   $$
   \alpha_n = 1 - \frac{p_n^{\mathrm{Cu|Au}}}{c_{\mathrm{Cu}}}
   \tag{1.3}
   $$

   where \(p_n^{\mathrm{Cu|Au}}\) is the probability that a site of shell \(n\) around an Au site is occupied by Cu, and \(c_{\mathrm{Cu}}\) is the overall Cu concentration (\(c_{\mathrm{Cu}} = 3/4\) in Cu₃Au).

3. Write the subprogram computing the long-range order parameter (function `ordre`).

   For the L1₂ structure of Cu₃Au, the long-range order parameter is defined by

   $$
   \eta = \frac{p_{\mathrm{Au}}^{(s)} - c_{\mathrm{Au}}}{1 - c_{\mathrm{Au}}}
   \tag{1.4}
   $$

   where \(p_{\mathrm{Au}}^{(s)}\) is the fraction of vertex sublattice sites occupied by Au atoms and \(c_{\mathrm{Au}} = 1/4\) is the overall Au concentration.

4. Test the complete program.

5. Determine the associated cohesive energy (at \(T = 0\ \mathrm{K}\)), starting from the ordered L1₂ structure.

6. Generate equilibrium configurations at different temperatures.

7. Calculate the different physical quantities as a function of temperature.

8. Determine the order-disorder transition temperature.

9. Discuss the results.

---

# Bibliographie

[1] Diu, Guttmann, Lederer, Roulet, *Physique statistique*, Hermann, complément IIIk p. 478, (1989).

[2] Y. Adda, J.M. Dupouy, J. Philibert, Y. Quéré, *Éléments de métallurgie physique*, INSTN-CEA Collection enseignement, vol. 3, 625–629 (1987).

[3] S.C. Moss, *J. Appl. Phys.* **35**, 3547 (1964).

[4] V. Rosato, M. Guillopé et B. Legrand, *Phil. Mag. A* **59**, 321 (1989).

[5] C. Rey-Losada, M. Hayoun, V. Pontikis, *Mat. Res. Soc. Symp. Proc.* **291**, 549 (1993).

---

*Hichem Dammak — Modélisation et Simulation Numérique*