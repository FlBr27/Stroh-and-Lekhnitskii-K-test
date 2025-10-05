# 6th-order Stroh and 4th-order Lekhnitskii K-tests for Grain Boundaries

This repository contains all necessary files to set up and run K-tests on generally oriented grain boundaries (GBs) in [LAMMPS](https://www.lammps.org). Thereby, two different variants are available: The 6th-order Stroh formalism should be used when the GB to be tested has no monoclinic material symmetry with respect to the out-of-plane direction. The 4th-order Lekhnitskii formalism can be used if the GB obeys the mentioned symmetry requirement. For the theory behind both approaches and further details, see  [Brunner et al. (2025) Modelling Simul. Mater. Sci. Eng. 33 035004](https://doi.org/10.1088/1361-651X/adba03). Please also cite this article when you are using the scripts provided in this repository.

The Stroh as well as the Lekhnitskii approach are used to compute displacement boundary conditions for the atomic simulation cell that correspond to a crack along a GB. To impose these boundary conditions in LAMMPS, we added two styles to the _displace_atoms_ LAMMPS command. The _crackaniso_ style implements the boundary conditions corresponding to the 4th-order Lekhnitskii formalism, while the _gbcrack_ style implements the 6th-order Stroh formalism. Two examples illustrating the use of the Stroh and Lekhnitskii formalism for K-tests of GBs are included in this repository.

Furthermore, a third example is provided, in which the results of a K-test are compared to the analytical predictions for brittle crack extensions (Griffith theory) and dislocation emission (Rice theory) within the framework of the STH formalism.  

# Installation

To enable the use of the _crackaniso_ and _gbcrack_ style of the _displace_atoms_ LAMMPS command, in the _src_ folder of your LAMMPS installation, replace the _displace_atoms.cpp_ file with the one contained in the _install_ folder. Afterwards, recompile your LAMMPS installation, and the new styles should be available.

# Examples 

The detailed description of the provided examples can be found in the corresponding sub-folder of the _examples_ folder. 

## A note on single crystals

Both the Lekhnitskii and Stroh formalism can be used for K-tests in single crystals as well. To this end, the crystal orientation is selected as being the same for the two interfacing grains. However, also in this case, care must be taken regarding the underlying symmetry requirements of both approaches. If a single crystal is to be simulated that does not exhibit monoclinic material symmetry with respect to the out-of-plane direction (see e.g. the single crystal example), the 6th-order Stroh formalism has to be used. Otherwise, also the 4th-order Lekhnitskii formalism can be used.      

# Concluding Remarks
The  _crackaniso_ style of the _displace_atoms_ LAMMPS command used in this repository is based on an earlier version by Predrag Andric (SKF) that was updated by Lei Zhang (Univerity of Groningen). Further details on the theory behind the 4th-order Lekhnitskii formalism and its application to fracture simulations are discussed in detail in: \
[Andric and Curtin (2019) Modelling Simul. Mater. Sci. Eng. 27, 013001](https://doi.org/10.1088/1361-651X/aae40c) \
Please also refer to this article when using the provided 4th-order Lekhnitskii approach.

If you find any bugs in the scripts of this repository or have suggestions for improvements or extensions, please feel free to contact me: \
Florian Brunner (f.brunner@rug.nl)
