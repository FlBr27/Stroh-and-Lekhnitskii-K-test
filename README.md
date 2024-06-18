# 6th-order Stroh and 4th-order Lekhnitskii K-tests for Grain Boundaries

This repository contains all necessary files to set up and run K-tests on generally oriented grain boundaries (GBs) in LAMMPS (https://www.lammps.org). Thereby, two different variants are available: The 6th-order Stroh formalism should be used when the GB to be tests has no monoclinic material symmetry with respect to the out-of-plane direction. The 4th-order Lekhnitskii formalism can be used if the GB obeys the mentioned symmetry requirement. For the theory behind both approaches and further details, see [link to the paper].

The Stroh as well as the Lekhnitskii approach are used to compute displacement boundary condition on the atomic simulation cell that correspond to a crack along a GB. To impose these boundary condition on the simulation cell in LAMMPS, we added two styles to the _displace_atoms_ LAMMPS command. The _crackaniso_ style implements the boundary conditions corresponding to the 4th-order Lekhnitskii formalism while the _gbcrack_ style implements the boundary conditions corresponding to the 6th-order Stroh formalism. Two example are included in this repository that illustrate the use of the Stroh and Lekhnitskii formalism for K-tests of GBs.   

# Installation

To enable the use of the _crackaniso_ and _gbcrack_ style of the _displace_atoms_ LAMMPS command, in the _src_ folder of your LAMMPS installation, replace the _displace_atoms.cpp_ file with the one contained in the "install" folder of this repository. Afterwards, recompile your LAMMPS installation and the new styles should be ready to use.

# Examples 

## Stroh K-test

The LAMMPS script _Stroh_K_Test.in_ in the _Stroh_example_ folder gives an example on how to run K-test simulations of GB configurations in LAMMPS, using the 6th-order Stroh formalism. In this script, a $\Sigma 3 \ (111) \ 109.5° \ [1\overline{1}0]$

