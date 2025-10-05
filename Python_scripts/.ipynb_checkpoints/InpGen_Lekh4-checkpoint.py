#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 16:48:40 2023

@author: Florian Brunner (f.brunner@rug.nl)

This script generates and exports the coefficients that are required as input
parameters for the "crackaniso" style of the "displace_atoms" command in LAMMPS
to generate the displacement field of a grain boundary according to the
4th-order Lekhnitskii (LSL) formalism.

"""

import numpy as np
import gb_funcs as gbfun

#################### Required Parameters ######################################

single_crystal = "no" # GB mode if "no" and single crystal mode if "yes"

# Grain boundary specification
sigma = 3                  # sigma value
phi_deg = 109.47122         # rotation angle in °
y1 = np.array([1,1,1])     # grain boundary normal w.r.t. grain 1
z1 = np.array([1,-1,0])      # common rotation axis

# Elastic Constants (Fe)
C_11 = 244.2857 # all in [GPa]
C_12 = 145.3217
C_44 = 116.3485

#################### Computation of Required Quantities #######################

if single_crystal == "no":

    # find the orientation vectors of both grains
    x1, x2, y2, z2 = gbfun.find_stgb_struc(phi_deg, y1, z1)
    
    # build the (unrotated) stiffness matrix
    C_el = gbfun.build_cubic_stiff(C_11, C_12, C_44)
    
    # find the stiffness matrices for grain 1 and 2
    C_el1 = gbfun.rotate_stiff_mat(C_el, x1, y1, z1)
    C_el2 = gbfun.rotate_stiff_mat(C_el, x2, y2, z2)
    
    # check if oscillatory behaviour is to be expected
    gbfun.Qu_Bassani_check(C_el1, C_el2)
    
    # Name of ouputfile (i.e. LAMMPS input file)
    outpf_name = 'GB_S'+str(sigma)+'_'+str(y1[0])+str(y1[1])+str(y1[2])+\
        '_'+str(z1[0])+str(z1[1])+str(z1[2])+'_Lekh4'
    
    coeff_Lekh4ps = gbfun.GB_InpGen_Lekh_4ps(C_el1, C_el2, write='LAMMPS',
                                        outpf_name=outpf_name,ex_hand='off')
    
elif single_crystal == "yes":
    
    x1 = np.cross(y1, z1)
    
    # build the (unrotated) stiffness matrix
    C_el = gbfun.build_cubic_stiff(C_11, C_12, C_44)
    
    # find the stiffness matrices for grain 1 and 2
    C_el1 = gbfun.rotate_stiff_mat(C_el, x1, y1, z1)
    
    # Name of ouputfile (i.e. LAMMPS input file)
    outpf_name = 'SC_'+str(y1[0])+str(y1[1])+str(y1[2])+\
        '_'+str(z1[0])+str(z1[1])+str(z1[2])+'_Lekh4'
    
    coeff_Lekh4ps = gbfun.GB_InpGen_Lekh_4ps(C_el1, C_el1, write='LAMMPS',
                                        outpf_name=outpf_name,ex_hand='off')
    
else:
    raise Exception("Single_crytal flag was incorrectly specified!")
