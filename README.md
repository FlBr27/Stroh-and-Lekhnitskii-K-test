# 6th-order Stroh and 4th-order Lekhnitskii K-tests for Grain Boundaries

This repository contains all necessary files to set up and run K-tests on generally oriented grain boundaries (GBs) in LAMMPS (https://www.lammps.org). Thereby, two different variants are available: The 6th-order Stroh formalism should be used when the GB to be tests has no monoclinic material symmetry with respect to the out-of-plane direction. The 4th-order Lekhnitskii formalism can be used if the GB obeys the mentioned symmetry requirement. For the theory behind both approaches and further details, see [link to the paper].

The Stroh as well as the Lekhnitskii approach are used to compute displacement boundary condition for the atomic simulation cell that correspond to a crack along a GB. To impose these boundary condition in LAMMPS, we added two styles to the _displace_atoms_ LAMMPS command. The _crackaniso_ style implements the boundary conditions corresponding to the 4th-order Lekhnitskii formalism while the _gbcrack_ style implements the 6th-order Stroh formalism. Two examples are included in this repository that illustrate the use of the Stroh and Lekhnitskii formalism for K-tests of GBs.   

# Installation

To enable the use of the _crackaniso_ and _gbcrack_ style of the _displace_atoms_ LAMMPS command, in the _src_ folder of your LAMMPS installation, replace the _displace_atoms.cpp_ file with the one contained in the _install_ folder. Afterwards, recompile your LAMMPS installation and the new styles should be available.

# Examples 

## 6th-order Stroh K-test

The LAMMPS script _Stroh_K_Test.in_ in the _Stroh_example_ folder illustrates how to run K-test simulations of GB configurations in LAMMPS, using the 6th-order Stroh formalism. In this script, a $\Sigma 9 \ (12\overline{2}) \ 152.7° [223]$ Fe GB is used as an example and loaded in mode I (pure tension) in $20$ steps with a stress intensity level ranging from $0.86$ to $1.26 \ \text{_MPa_}\sqrt{m}$. Thereby, the _Fe_mm.eam.fs_ interatomic potential is used for Fe, since it is by default included in the _potentials_ folder of LAMMPS. It also has to be noted that this GB has no monoclinic material symmetry with respect to the out-of-plane direction and therefore really requires the use of the 6th-order Stroh formalism. 

The script reads the (fully relaxed) simulation cell with the GB from the _corr_optimum_gb_S9.restart_ restart file. The simulation cell can be examined e.g. with OVITO (https://www.ovito.org/) by opening the _corr_optimum_gb_S9.data_ file.

In line $146$ and $158$ the _displace_atoms_ command is called using the _gbcrack_ style. This style requires the specification of $67$ input parameters (mostly real and imaginary parts of matrices and vectors). However, in general, only the stress intensity factors $K_{I}$, $K_{II}$, and $K_{III}$ as well as the coordinates of the (initial) crack tip position $x_{tip}$ and $y_{tip}$ have to be provided by the user. The remaining parameters can be calculated with the Python script _InpGen_Stroh6.py_ in the _Python_scripts_ folder and imported in LAMMPS using the _include_ command. This is shown in line $103$ in _Stroh_K_Test.in_ with the _GB_S9_12-2_223_Stroh6_ parameter file. More information on the input parameters of the _gbcrack_ style can be found in the _how_to_gbcrack.txt_ file, also contained in the _Stroh_example_ folder.

Note: To simulate a different GB (material, character, orientation, ...) a new simulation cell restart file and a new input parameter file are required. 

When executed with LAMMPS, _Stroh_K_Test.in_ produces an output (dump) file for each loading step that can be inspected e.g. in OVITO. If the script was run correctly, one should see brittle crack extension of the $\Sigma 9 \ (12\overline{2}) \ 152.7° [223]$ GB with a first larger crack length increase from loading step $9$ to $10$. 

## 4th-order Lekhnitskii K-test

The LAMMPS script _Lekh_K_Test.in_ in the _Lekhnitskii_example_ folder illustrates how to run K-test simulations of GB configurations in LAMMPS, using the 4th-order Lekhnitskii formalism. In this script, a $\Sigma 3 \ (111) \ 109.5° [1\overline{1}0]$ Fe GB is used as an example and loaded in mode I (pure tension) in $20$ steps with a stress intensity level ranging from $0.80$ to $1.20 MPa\sqrt{m}$. Thereby, the _Fe_mm.eam.fs_ interatomic potential is again used. It also has to be noted that this GB has monoclinic material symmetry with respect to the out-of-plane direction and therefore the 4th-order Lekhnitskii formalism can be used. 

The script reads the (fully relaxed) simulation cell with the GB from the _corr_optimum_gb_S3.restart_ restart file. The simulation cell can again be examined e.g. with OVITO (https://www.ovito.org/) by opening the _corr_optimum_gb_S3.data_ file.

In line $147$ and $151$ the _displace_atoms_ command is called using the _crackaniso_ style. This style requires the specification of $15$ input parameters (mostly real and imaginary parts of scalars). However, in general, only the stress intensity factor $K_{I}$ and the coordinates of the (initial) crack tip position $x_{tip}$ and $y_{tip}$ have to be provided by the user, since the _crackaniso_ style is limited to mode I loading only. The remaining parameters can be calculated with the Python script _InpGen_Lekh4.py_ in the _Python_scripts_ folder and imported in LAMMPS using the _include_ command. This is shown in line $104$ in _Lekh_K_Test.in_ with the _GB_S3_111_1-10_Lekh4_ parameter file. More information regarding the input parameters of the _crackaniso_ style can be found in the _how_to_crackaniso.txt_ file, also contained in the _Lekhnitskii_example_ folder.

Note 1: To simulate a different GB (material, character, orientation, ...) a new simulation cell restart file and a new input parameter file are required. 
Note 2: When providing stress intensity values for the _crackaniso_ style, they have to be specified in units of $[100 MPa\sqrt{m}]$ contrary to the _gbcrack_ style where stress intensities are specified in units of $[MPa\sqrt{m}]$. 

When executed with LAMMPS, _Lekh_K_Test.in_ produces an output (dump) file for each loading step that can be inspected e.g. in OVITO. If the script was run correctly, one should see twinning fault emission from the crack tip that starts at load step $16$.

## A note on single crystals

Both the Lekhnitskii and Stroh formalism can be used for K-tests in single crystals as well. To this end, the crystal orientation is simply selected as being the same for the two interfacing grains. However, also in this case, care must be taken regarding the underlying symmetry requirements of both approaches: If a single crystal is to be simulated that does not exhibit monoclinic material symmetry with respect to the out-of-plane direction, the 6th-order Stroh formalism has to be used. Otherwise, the 4th-order Lekhnitskii formalism can be used too.      

# Concluding Remarks
If you find any bugs in the scripts of this repository or have suggestions for improvements or extensions, please feel free to contact me:
Florian Brunner (f.brunner@rug.nl)
