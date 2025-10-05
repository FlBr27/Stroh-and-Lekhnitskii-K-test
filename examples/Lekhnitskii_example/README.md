## 4th-order Lekhnitskii K-test

The LAMMPS script _Lekh_K_Test.in_ illustrates how to run K-test simulations of GB configurations in LAMMPS, using the 4th-order Lekhnitskii formalism. 
In this script, a $\Sigma 3 \ (111) \ 109.5° [1\overline{1}0]$ Fe GB is used as an example and loaded in mode I (pure tension) in $20$ steps with a stress intensity 
level ranging from $0.80$ to $1.20 \ \text{MPa}\sqrt{m}$. Thereby, the _Fe_mm.eam.fs_ interatomic potential is used for Fe, since it is by default included in the _potentials_ 
folder of LAMMPS. It also has to be noted that this GB has monoclinic material symmetry with respect to the out-of-plane direction, and therefore, the 4th-order Lekhnitskii 
formalism can be used.

The script reads the (fully relaxed) simulation cell with the GB from the _corr_optimum_gb_S3.restart_ restart file. The simulation cell can be examined, e.g. with 
[OVITO](https://www.ovito.org/) by opening the _corr_optimum_gb_S3.data_ file.

In lines $147$ and $151$, the _displace_atoms_ command is called using the _crackaniso_ style. This style requires the specification of $15$ input parameters (mostly real and 
imaginary parts of scalars). However, in general, only the stress intensity factor $K_{I}$ and the coordinates of the (initial) crack tip position $x_{tip}$ and $y_{tip}$ have 
to be provided by the user, since the _crackaniso_ style is limited to mode I loading only. The remaining parameters can be calculated with the Jupyter Notebook script 
_InpGen_Lekh4.ipynb_, provided in the _python_scripts_ folder and imported in LAMMPS using the _include_ command. This is shown in line $104$ in _Lekh_K_Test.in_ with 
the _GB_S3_111_1-10_Lekh4_ parameter file. More information regarding the input parameters of the _crackaniso_ style can be found in the _how_to_crackaniso.txt_ file.

**Note 1:** To simulate a different GB (material, character, orientation, ...), a new simulation cell restart file and a new input parameter file are required. \
**Note 2:** When providing stress intensity values for the _crackaniso_ style, they have to be specified in units of $[100 \ \text{MPa}\sqrt{m}]$, contrary to the _gbcrack_ 
style where stress intensities are specified in units of $[\text{MPa}\sqrt{m}]$.

When executed with LAMMPS, _Lekh_K_Test.in_ produces an output (dump) file for each loading step that can be inspected, e.g. in OVITO. If the script was run correctly, one should 
see twinning fault emission from the crack tip that starts at load step $16$.
