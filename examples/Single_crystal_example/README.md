## Comparison of a single crystal K-test with theory

The purpose of this example is to demonstrate how the predictions of a K-test simulation can be compared to analytical theories. Two of the most important theories in this regard are, on the one hand, the theory of brittle crack extension, devised by Griffith (Phil. Trans. Royal Soc. A 221, pp. 163–198 (1920)). On the other hand, 
Rice (J. Mech. Phys. Solids Vol. 40, No. 2, pp. 239-271 (1992)) proposed a theory for dislocation emission from a crack tip that was subsequently extended to anisotropic cases
by Sun & Beltz (J. Mech. Phys. Solids Vol. 42, No. 12, pp. 1905-1932 (1994)). These two theories can be considered to validate the predictions of K-tests, as exemplified hereafter.

The LAMMPS script _Stroh_K_Test.in_ is used to perform a K-test simulation of a single crystal configuration in LAMMPS, using the 6th-order Stroh formalism. 
In this script, a $(111) \ [11\overline{2}]$ Fe single crystal is considered and loaded in mode I (pure tension) in $20$ steps with a stress intensity level 
ranging from $0.95$ to $1.35 \ \text{MPa}\sqrt{m}$. Thereby, the _Fe_mm.eam.fs_ interatomic potential is used for Fe, since it is by default included in the _potentials_ 
folder of LAMMPS. It also has to be noted that this single crystal has no monoclinic material symmetry with respect to the out-of-plane direction and therefore requires the use 
of the 6th-order Stroh formalism. Further details on the K-test setup can be found in the _Lekhnitskii example_ and the _Stroh example_. 
 
When executed with LAMMPS, _Stroh_K_Test.in_ produces an output (dump) file for each loading step, which can be inspected, e.g., in OVITO. If the script was run correctly, 
one should see the emission of an edge dislocation at an angle $\theta = 90°$ from the crack plane at load step $11$, i.e. for a stress intensity level of $1.17 \ \text{MPa}\sqrt{m}$. 

To assess whether this prediction is reasonable or not, the Griffith value $K_{G}$ for brittle crack extension is computed using the Jupyter Notebook script _Griffith_compute.ipynb_ in the _python_scripts_ folder. Executing this script for the here treated case yields $K_{G} = 0.981\ \text{MPa}\sqrt{m}$. On the other hand, 
the relevant slip systems for dislocation emission and their corresponding Rice value $K_{Ie}$ are determined via the _Rice_compute.ipynb_ Jupyter Notebook script, 
also to be found in the  _python_scripts_ folder. By executing the script, one finds that the most favourable slip systems in this case are the $\\{110\\}\langle111\rangle$ systems, which are
inclined to the crack plane (x-z-plane) by an angle $\theta = \pm 90°$, and to the crack front plane (x-y-plane) by an angle $\phi = 0°$.

Comparing the theoretical predictions, brittle fracture would be the expected outcome of the conducted K-test since $K_{G} < K_{Ie}$. However, dislocation emission is observed
in the slip system predicted by the Rice analysis. This discrepancy is most likely the result of the used interatomic potential, which is a (simple)
EAM model. This issue was discussed in detail by e.g. Möller & Bitzek (Modelling Simul. Mater. Sci. Eng. 22 (2014) 045002).



