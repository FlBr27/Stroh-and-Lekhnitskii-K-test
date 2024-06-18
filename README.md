# 6th-order Stroh and 4th-order Lekhnitskii K-tests for Grain Boundaries

This repository contains all necessary files to set up and run K-tests on generally oriented grain boundaries (GBs) in LAMMPS (https://www.lammps.org). Thereby, two different variants are available: The 6th-order Stroh formalism should be used when the GB to be tests has no monoclinic material symmetry with respect to the out-of-plane direction. The 4th-order Lekhnitskii formalism can be used if the GB obeys the mentioned symmetry requirement. For the theory behind both approaches and further details, see [link to the paper].

The Stroh as well as the Lekhnitskii approach are used to compute displacement boundary condition on the atomic simulation cell that correspond to a crack along a GB. To impose these boundary condition on the simulation cell in LAMMPS, we added two styles to the "displace_atoms" LAMMPS command. The "crackaniso" style implements the boundary conditions corresponding to the 4th-order Lekhnitskii formalism while the "gbcrack" style implements the boundary conditions corresponding to the 6th-order Stroh formalism. Two example are included in this repository that illustrate the use of the Stroh and Lekhnitskii formalism for K-tests of GBs.   

# Installation

To enable the use of the "crackaniso" and "gbcrack" style of the "displace_atoms" LAMMPS command, in the "src" folder of your LAMMPS installation, replace the "displace_atoms.cpp" file with the one contained in the "install" folder of this repository. Afterwards, recompile your LAMMPS installation and the new styles should be ready to use.

# Examples 

## Stroh K-test

In the 
