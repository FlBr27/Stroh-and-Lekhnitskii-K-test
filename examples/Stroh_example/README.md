## 6th-order Stroh K-test

The LAMMPS script _Stroh_K_Test.in_ illustrates how to run K-test simulations of GB configurations in LAMMPS, using the 6th-order Stroh formalism. 
In this script, a $\Sigma 9 \ (12\overline{2}) \ 152.7° [223]$ Fe GB is used as an example and loaded in mode I (pure tension) in $20$ steps with a stress intensity level 
ranging from $0.86$ to $1.26 \ \text{MPa}\sqrt{m}$. Thereby, the _Fe_mm.eam.fs_ interatomic potential is used for Fe, since it is by default included in the _potentials_ folder 
of LAMMPS. It also has to be noted that this GB has no monoclinic material symmetry with respect to the out-of-plane direction and therefore requires the use of the 
6th-order Stroh formalism. 

The script reads the (fully relaxed) simulation cell with the GB from the _corr_optimum_gb_S9.restart_ restart file. The simulation cell can be examined, e.g. with 
[OVITO](https://www.ovito.org/) by opening the _corr_optimum_gb_S9.data_ file.

In lines $146$ and $158$, the _displace_atoms_ command is called using the _gbcrack_ style. This style requires the specification of $67$ input parameters (mostly real and 
imaginary parts of matrices and vectors). However, in general, only the stress intensity factors $K_{I}$, $K_{II}$, and $K_{III}$ as well as the coordinates of the (initial) 
crack tip position $x_{tip}$ and $y_{tip}$ have to be provided by the user. The remaining parameters can be calculated with the Jupyter Notebook script _InpGen_Stroh6.ipynb_ 
in the _python_scripts_ folder, and imported in LAMMPS using the _include_ command. This is shown in line $103$ in _Stroh_K_Test.in_ with the _GB_S9_12-2_223_Stroh6_ parameter file. More information on the input parameters of the _gbcrack_ style can be found in the _how_to_gbcrack.txt_ file.

**Note:** To simulate a different GB (material, character, orientation, ...), a new simulation cell restart file and a new input parameter file are required.

When executed with LAMMPS, _Stroh_K_Test.in_ produces an output (dump) file for each loading step, which can be inspected, e.g. in OVITO. If the script was run correctly, 
one should see brittle crack extension of the $\Sigma 9 \ (12\overline{2}) \ 152.7° [223]$ GB with a first larger crack length increase from loading step $9$ to $10$. 
