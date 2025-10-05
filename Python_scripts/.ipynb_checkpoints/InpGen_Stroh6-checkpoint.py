# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 20:25:59 2023

@author: Florian Brunner (f.brunner@rug.nl)

This script generates and exports the coefficients that are required as input
parameters for the "gbcrack" style of the "displace_atoms" command in LAMMPS
to generate the displacement field of a grain boundary according to the
6th-order Stroh (STH) formalism.
"""

import os
import numpy as np
from numpy import linalg as la
import gb_funcs as gbfun

#################### Required Parameters ######################################

single_crystal = "no" # GB mode if "no" and single crystal mode if "yes"

# Grain boundary specification
sigma = 9                  # sigma value
phi_deg = 152.7340         # rotation angle in °
y1 = np.array([1,2,-2])     # grain boundary normal w.r.t. grain 1
z1 = np.array([2,2,3])      # common rotation axis

# Elastic Constants (Fe)
C_11 = 244.2857 # all in [GPa]
C_12 = 145.3217
C_44 = 116.3485

# Crack length and grain (region) to be simulated
a = 10000.0 # [A]

#################### Computation of Required Quantities #######################

if single_crystal == "no":

    # Name of ouputfile (i.e. LAMMPS input file)
    outpf_name = 'GB_S'+str(sigma)+'_'+str(y1[0])+str(y1[1])+str(y1[2])+\
        '_'+str(z1[0])+str(z1[1])+str(z1[2])+'_Stroh6'
    
    # find the orientation vectors of both grains
    x1, x2, y2, z2 = gbfun.find_stgb_struc(phi_deg, y1, z1)
    
    # build the (unrotated) stiffness matrix
    C_el = gbfun.build_cubic_stiff(C_11, C_12, C_44)
    
    # find the stiffness matrices for grain 1 and 2
    C_el1 = gbfun.rotate_stiff_mat(C_el, x1, y1, z1)
    C_el2 = gbfun.rotate_stiff_mat(C_el, x2, y2, z2)
    
    # calculate p, A and B of the Stroh formalism
    p1, A1, B1, = gbfun.Stroh6(C_el1)
    p2, A2, B2, = gbfun.Stroh6(C_el2)
    
    B1_inv = la.inv(B1)
    B2_inv = la.inv(B2)
    
    # calculate the required parameters of the Stroh formalism
    Param = gbfun.Stroh6_Coeff(p1, A1, B1, p2, A2, B2, a)
    
    S_brev = Param['S_brev']
    Phi_inv = Param['Phi_inv']
    
    # determing E_11 for the Griffith criterion
    # inverse impedance tensors
    M_inv1 = 1j*A1 @ B1_inv
    M_inv2 = 1j*A2 @ B2_inv
    
    # bimaterial matrix
    M_st = M_inv1 + np.conjugate(M_inv2)
    
    D = np.real(M_st)
    W = -np.imag(M_st)
    
    E = D + W @ la.inv(D) @ W
    
    # write coefficients to LAMMPS readable file
    
    # writing parameters to file
    outpf_name = './outputs/' + outpf_name
    os.makedirs(os.path.dirname(outpf_name), exist_ok=True)
    f = open(outpf_name, 'w')
    
    f.write('# Input parameters for the 6th order Stroh K-test \n')
    f.write('\n')
    
    # crack length a
    f.write('variable a equal {:.3f} \n'.format(a))
    f.write('\n')
    
    f.write('### Parameters region 1\n')
    
    # p_j's
    f.write('# p1_j\n')
    f.write('variable p11_re equal {:.9f} \n'.format(p1[0].real))
    f.write('variable p11_im equal {:.9f} \n'.format(p1[0].imag))
    
    f.write('variable p12_re equal {:.9f} \n'.format(p1[1].real))
    f.write('variable p12_im equal {:.9f} \n'.format(p1[1].imag))
    
    f.write('variable p13_re equal {:.9f} \n'.format(p1[2].real))
    f.write('variable p13_im equal {:.9f} \n'.format(p1[2].imag))
    
    # A-matrix
    f.write('# A1-matrix\n')
    f.write('variable a111_re equal {:.9f} \n'.format(A1[0,0].real))
    f.write('variable a111_im equal {:.9f} \n'.format(A1[0,0].imag))
    
    f.write('variable a112_re equal {:.9f} \n'.format(A1[0,1].real))
    f.write('variable a112_im equal {:.9f} \n'.format(A1[0,1].imag))
    
    f.write('variable a113_re equal {:.9f} \n'.format(A1[0,2].real))
    f.write('variable a113_im equal {:.9f} \n'.format(A1[0,2].imag))
    
    f.write('variable a121_re equal {:.9f} \n'.format(A1[1,0].real))
    f.write('variable a121_im equal {:.9f} \n'.format(A1[1,0].imag))
    
    f.write('variable a122_re equal {:.9f} \n'.format(A1[1,1].real))
    f.write('variable a122_im equal {:.9f} \n'.format(A1[1,1].imag))
    
    f.write('variable a123_re equal {:.9f} \n'.format(A1[1,2].real))
    f.write('variable a123_im equal {:.9f} \n'.format(A1[1,2].imag))
    
    f.write('variable a131_re equal {:.9f} \n'.format(A1[2,0].real))
    f.write('variable a131_im equal {:.9f} \n'.format(A1[2,0].imag))
    
    f.write('variable a132_re equal {:.9f} \n'.format(A1[2,1].real))
    f.write('variable a132_im equal {:.9f} \n'.format(A1[2,1].imag))
    
    f.write('variable a133_re equal {:.9f} \n'.format(A1[2,2].real))
    f.write('variable a133_im equal {:.9f} \n'.format(A1[2,2].imag))
    
    # B_inv-matrix
    f.write('# B1_inv-matrix\n')
    f.write('variable b111_re equal {:.9f} \n'.format(B1_inv[0,0].real))
    f.write('variable b111_im equal {:.9f} \n'.format(B1_inv[0,0].imag))
    
    f.write('variable b112_re equal {:.9f} \n'.format(B1_inv[0,1].real))
    f.write('variable b112_im equal {:.9f} \n'.format(B1_inv[0,1].imag))
    
    f.write('variable b113_re equal {:.9f} \n'.format(B1_inv[0,2].real))
    f.write('variable b113_im equal {:.9f} \n'.format(B1_inv[0,2].imag))
    
    f.write('variable b121_re equal {:.9f} \n'.format(B1_inv[1,0].real))
    f.write('variable b121_im equal {:.9f} \n'.format(B1_inv[1,0].imag))
    
    f.write('variable b122_re equal {:.9f} \n'.format(B1_inv[1,1].real))
    f.write('variable b122_im equal {:.9f} \n'.format(B1_inv[1,1].imag))
    
    f.write('variable b123_re equal {:.9f} \n'.format(B1_inv[1,2].real))
    f.write('variable b123_im equal {:.9f} \n'.format(B1_inv[1,2].imag))
    
    f.write('variable b131_re equal {:.9f} \n'.format(B1_inv[2,0].real))
    f.write('variable b131_im equal {:.9f} \n'.format(B1_inv[2,0].imag))
    
    f.write('variable b132_re equal {:.9f} \n'.format(B1_inv[2,1].real))
    f.write('variable b132_im equal {:.9f} \n'.format(B1_inv[2,1].imag))
    
    f.write('variable b133_re equal {:.9f} \n'.format(B1_inv[2,2].real))
    f.write('variable b133_im equal {:.9f} \n'.format(B1_inv[2,2].imag))
    
    f.write('\n')
    f.write('### Parameters region 2\n')
    
    # p_j's
    f.write('# p2_j\n')
    f.write('variable p21_re equal {:.9f} \n'.format(p2[0].real))
    f.write('variable p21_im equal {:.9f} \n'.format(p2[0].imag))
    
    f.write('variable p22_re equal {:.9f} \n'.format(p2[1].real))
    f.write('variable p22_im equal {:.9f} \n'.format(p2[1].imag))
    
    f.write('variable p23_re equal {:.9f} \n'.format(p2[2].real))
    f.write('variable p23_im equal {:.9f} \n'.format(p2[2].imag))
    
    # A-matrix
    f.write('# A2-matrix\n')
    f.write('variable a211_re equal {:.9f} \n'.format(A2[0,0].real))
    f.write('variable a211_im equal {:.9f} \n'.format(A2[0,0].imag))
    
    f.write('variable a212_re equal {:.9f} \n'.format(A2[0,1].real))
    f.write('variable a212_im equal {:.9f} \n'.format(A2[0,1].imag))
    
    f.write('variable a213_re equal {:.9f} \n'.format(A2[0,2].real))
    f.write('variable a213_im equal {:.9f} \n'.format(A2[0,2].imag))
    
    f.write('variable a221_re equal {:.9f} \n'.format(A2[1,0].real))
    f.write('variable a221_im equal {:.9f} \n'.format(A2[1,0].imag))
    
    f.write('variable a222_re equal {:.9f} \n'.format(A2[1,1].real))
    f.write('variable a222_im equal {:.9f} \n'.format(A2[1,1].imag))
    
    f.write('variable a223_re equal {:.9f} \n'.format(A2[1,2].real))
    f.write('variable a223_im equal {:.9f} \n'.format(A2[1,2].imag))
    
    f.write('variable a231_re equal {:.9f} \n'.format(A2[2,0].real))
    f.write('variable a231_im equal {:.9f} \n'.format(A2[2,0].imag))
    
    f.write('variable a232_re equal {:.9f} \n'.format(A2[2,1].real))
    f.write('variable a232_im equal {:.9f} \n'.format(A2[2,1].imag))
    
    f.write('variable a233_re equal {:.9f} \n'.format(A2[2,2].real))
    f.write('variable a233_im equal {:.9f} \n'.format(A2[2,2].imag))
    
    # B_inv-matrix
    f.write('# B2_inv-matrix\n')
    f.write('variable b211_re equal {:.9f} \n'.format(B2_inv[0,0].real))
    f.write('variable b211_im equal {:.9f} \n'.format(B2_inv[0,0].imag))
    f.write('variable b212_re equal {:.9f} \n'.format(B2_inv[0,1].real))
    f.write('variable b212_im equal {:.9f} \n'.format(B2_inv[0,1].imag))
    f.write('variable b213_re equal {:.9f} \n'.format(B2_inv[0,2].real))
    f.write('variable b213_im equal {:.9f} \n'.format(B2_inv[0,2].imag))
    f.write('variable b221_re equal {:.9f} \n'.format(B2_inv[1,0].real))
    f.write('variable b221_im equal {:.9f} \n'.format(B2_inv[1,0].imag))
    f.write('variable b222_re equal {:.9f} \n'.format(B2_inv[1,1].real))
    f.write('variable b222_im equal {:.9f} \n'.format(B2_inv[1,1].imag))
    f.write('variable b223_re equal {:.9f} \n'.format(B2_inv[1,2].real))
    f.write('variable b223_im equal {:.9f} \n'.format(B2_inv[1,2].imag))
    f.write('variable b231_re equal {:.9f} \n'.format(B2_inv[2,0].real))
    f.write('variable b231_im equal {:.9f} \n'.format(B2_inv[2,0].imag))
    f.write('variable b232_re equal {:.9f} \n'.format(B2_inv[2,1].real))
    f.write('variable b232_im equal {:.9f} \n'.format(B2_inv[2,1].imag))
    f.write('variable b233_re equal {:.9f} \n'.format(B2_inv[2,2].real))
    f.write('variable b233_im equal {:.9f} \n'.format(B2_inv[2,2].imag))
    
    f.write('\n')
    f.write('# Interfacial parameters\n')
    
    # S_breve-matrix
    f.write('# S_breve-matrix\n')
    
    f.write('variable s_11 equal {:.9f} \n'.format(S_brev[0,0]))
    f.write('variable s_12 equal {:.9f} \n'.format(S_brev[0,1]))
    f.write('variable s_13 equal {:.9f} \n'.format(S_brev[0,2]))
    f.write('variable s_21 equal {:.9f} \n'.format(S_brev[1,0]))
    f.write('variable s_22 equal {:.9f} \n'.format(S_brev[1,1]))
    f.write('variable s_23 equal {:.9f} \n'.format(S_brev[1,2]))
    f.write('variable s_31 equal {:.9f} \n'.format(S_brev[2,0]))
    f.write('variable s_32 equal {:.9f} \n'.format(S_brev[2,1]))
    f.write('variable s_33 equal {:.9f} \n'.format(S_brev[2,2]))
    
    # Phi_inv-matrix
    f.write('# Phi_inv-matrix\n')
    f.write('variable Phi_11 equal {:.9f} \n'.format(Phi_inv[0,0]))
    f.write('variable Phi_12 equal {:.9f} \n'.format(Phi_inv[0,1]))
    f.write('variable Phi_13 equal {:.9f} \n'.format(Phi_inv[0,2]))
    f.write('variable Phi_21 equal {:.9f} \n'.format(Phi_inv[1,0]))
    f.write('variable Phi_22 equal {:.9f} \n'.format(Phi_inv[1,1]))
    f.write('variable Phi_23 equal {:.9f} \n'.format(Phi_inv[1,2]))
    f.write('variable Phi_31 equal {:.9f} \n'.format(Phi_inv[2,0]))
    f.write('variable Phi_32 equal {:.9f} \n'.format(Phi_inv[2,1]))
    f.write('variable Phi_33 equal {:.9f} \n'.format(Phi_inv[2,2]))
    
    f.write('\n')
    f.write('variable E_22 equal {:.9f} \n'.format(E[1,1]))
    
elif single_crystal == "yes":

    # Name of ouputfile (i.e. LAMMPS input file)
    outpf_name = 'SC_'+str(y1[0])+str(y1[1])+str(y1[2])+\
        '_'+str(z1[0])+str(z1[1])+str(z1[2])+'_Stroh6'
    
    x1 = np.cross(y1, z1)
    
    # build the (unrotated) stiffness matrix
    C_el = gbfun.build_cubic_stiff(C_11, C_12, C_44)
    
    # find the stiffness matrices for grain 1 and 2
    C_el1 = gbfun.rotate_stiff_mat(C_el, x1, y1, z1)
    
    # calculate p, A and B of the Stroh formalism
    p1, A1, B1, = gbfun.Stroh6(C_el1)
    
    B1_inv = la.inv(B1)
    
    # calculate the required parameters of the Stroh formalism
    Param = gbfun.Stroh6_Coeff(p1, A1, B1, p1, A1, B1, a)
    
    S_brev = Param['S_brev']
    Phi_inv = Param['Phi_inv']
    
    # determing E_11 for the Griffith criterion
    # inverse impedance tensors
    M_inv1 = 1j*A1 @ B1_inv
    
    # bimaterial matrix
    M_st = M_inv1 + np.conjugate(M_inv1)
    
    D = np.real(M_st)
    W = -np.imag(M_st)
    
    E = D + W @ la.inv(D) @ W
    
    # write coefficients to LAMMPS readable file
    
    # writing parameters to file
    outpf_name = './outputs/' + outpf_name
    os.makedirs(os.path.dirname(outpf_name), exist_ok=True)
    f = open(outpf_name, 'w')
    
    f.write('# Input parameters for the 6th order Stroh K-test \n')
    f.write('\n')
    
    # crack length a
    f.write('variable a equal {:.3f} \n'.format(a))
    f.write('\n')
    
    f.write('### Parameters region 1\n')
    
    # p_j's
    f.write('# p1_j\n')
    f.write('variable p11_re equal {:.9f} \n'.format(p1[0].real))
    f.write('variable p11_im equal {:.9f} \n'.format(p1[0].imag))
    
    f.write('variable p12_re equal {:.9f} \n'.format(p1[1].real))
    f.write('variable p12_im equal {:.9f} \n'.format(p1[1].imag))
    
    f.write('variable p13_re equal {:.9f} \n'.format(p1[2].real))
    f.write('variable p13_im equal {:.9f} \n'.format(p1[2].imag))
    
    # A-matrix
    f.write('# A1-matrix\n')
    f.write('variable a111_re equal {:.9f} \n'.format(A1[0,0].real))
    f.write('variable a111_im equal {:.9f} \n'.format(A1[0,0].imag))
    
    f.write('variable a112_re equal {:.9f} \n'.format(A1[0,1].real))
    f.write('variable a112_im equal {:.9f} \n'.format(A1[0,1].imag))
    
    f.write('variable a113_re equal {:.9f} \n'.format(A1[0,2].real))
    f.write('variable a113_im equal {:.9f} \n'.format(A1[0,2].imag))
    
    f.write('variable a121_re equal {:.9f} \n'.format(A1[1,0].real))
    f.write('variable a121_im equal {:.9f} \n'.format(A1[1,0].imag))
    
    f.write('variable a122_re equal {:.9f} \n'.format(A1[1,1].real))
    f.write('variable a122_im equal {:.9f} \n'.format(A1[1,1].imag))
    
    f.write('variable a123_re equal {:.9f} \n'.format(A1[1,2].real))
    f.write('variable a123_im equal {:.9f} \n'.format(A1[1,2].imag))
    
    f.write('variable a131_re equal {:.9f} \n'.format(A1[2,0].real))
    f.write('variable a131_im equal {:.9f} \n'.format(A1[2,0].imag))
    
    f.write('variable a132_re equal {:.9f} \n'.format(A1[2,1].real))
    f.write('variable a132_im equal {:.9f} \n'.format(A1[2,1].imag))
    
    f.write('variable a133_re equal {:.9f} \n'.format(A1[2,2].real))
    f.write('variable a133_im equal {:.9f} \n'.format(A1[2,2].imag))
    
    # B_inv-matrix
    f.write('# B1_inv-matrix\n')
    f.write('variable b111_re equal {:.9f} \n'.format(B1_inv[0,0].real))
    f.write('variable b111_im equal {:.9f} \n'.format(B1_inv[0,0].imag))
    
    f.write('variable b112_re equal {:.9f} \n'.format(B1_inv[0,1].real))
    f.write('variable b112_im equal {:.9f} \n'.format(B1_inv[0,1].imag))
    
    f.write('variable b113_re equal {:.9f} \n'.format(B1_inv[0,2].real))
    f.write('variable b113_im equal {:.9f} \n'.format(B1_inv[0,2].imag))
    
    f.write('variable b121_re equal {:.9f} \n'.format(B1_inv[1,0].real))
    f.write('variable b121_im equal {:.9f} \n'.format(B1_inv[1,0].imag))
    
    f.write('variable b122_re equal {:.9f} \n'.format(B1_inv[1,1].real))
    f.write('variable b122_im equal {:.9f} \n'.format(B1_inv[1,1].imag))
    
    f.write('variable b123_re equal {:.9f} \n'.format(B1_inv[1,2].real))
    f.write('variable b123_im equal {:.9f} \n'.format(B1_inv[1,2].imag))
    
    f.write('variable b131_re equal {:.9f} \n'.format(B1_inv[2,0].real))
    f.write('variable b131_im equal {:.9f} \n'.format(B1_inv[2,0].imag))
    
    f.write('variable b132_re equal {:.9f} \n'.format(B1_inv[2,1].real))
    f.write('variable b132_im equal {:.9f} \n'.format(B1_inv[2,1].imag))
    
    f.write('variable b133_re equal {:.9f} \n'.format(B1_inv[2,2].real))
    f.write('variable b133_im equal {:.9f} \n'.format(B1_inv[2,2].imag))
    
    f.write('\n')
    f.write('### Parameters region 2\n')
    
    # p_j's
    f.write('# p2_j\n')
    f.write('variable p21_re equal {:.9f} \n'.format(p1[0].real))
    f.write('variable p21_im equal {:.9f} \n'.format(p1[0].imag))
    
    f.write('variable p22_re equal {:.9f} \n'.format(p1[1].real))
    f.write('variable p22_im equal {:.9f} \n'.format(p1[1].imag))
    
    f.write('variable p23_re equal {:.9f} \n'.format(p1[2].real))
    f.write('variable p23_im equal {:.9f} \n'.format(p1[2].imag))
    
    # A-matrix
    f.write('# A2-matrix\n')
    f.write('variable a211_re equal {:.9f} \n'.format(A1[0,0].real))
    f.write('variable a211_im equal {:.9f} \n'.format(A1[0,0].imag))
    
    f.write('variable a212_re equal {:.9f} \n'.format(A1[0,1].real))
    f.write('variable a212_im equal {:.9f} \n'.format(A1[0,1].imag))
    
    f.write('variable a213_re equal {:.9f} \n'.format(A1[0,2].real))
    f.write('variable a213_im equal {:.9f} \n'.format(A1[0,2].imag))
    
    f.write('variable a221_re equal {:.9f} \n'.format(A1[1,0].real))
    f.write('variable a221_im equal {:.9f} \n'.format(A1[1,0].imag))
    
    f.write('variable a222_re equal {:.9f} \n'.format(A1[1,1].real))
    f.write('variable a222_im equal {:.9f} \n'.format(A1[1,1].imag))
    
    f.write('variable a223_re equal {:.9f} \n'.format(A1[1,2].real))
    f.write('variable a223_im equal {:.9f} \n'.format(A1[1,2].imag))
    
    f.write('variable a231_re equal {:.9f} \n'.format(A1[2,0].real))
    f.write('variable a231_im equal {:.9f} \n'.format(A1[2,0].imag))
    
    f.write('variable a232_re equal {:.9f} \n'.format(A1[2,1].real))
    f.write('variable a232_im equal {:.9f} \n'.format(A1[2,1].imag))
    
    f.write('variable a233_re equal {:.9f} \n'.format(A1[2,2].real))
    f.write('variable a233_im equal {:.9f} \n'.format(A1[2,2].imag))
    
    # B_inv-matrix
    f.write('# B2_inv-matrix\n')
    f.write('variable b211_re equal {:.9f} \n'.format(B1_inv[0,0].real))
    f.write('variable b211_im equal {:.9f} \n'.format(B1_inv[0,0].imag))
    f.write('variable b212_re equal {:.9f} \n'.format(B1_inv[0,1].real))
    f.write('variable b212_im equal {:.9f} \n'.format(B1_inv[0,1].imag))
    f.write('variable b213_re equal {:.9f} \n'.format(B1_inv[0,2].real))
    f.write('variable b213_im equal {:.9f} \n'.format(B1_inv[0,2].imag))
    f.write('variable b221_re equal {:.9f} \n'.format(B1_inv[1,0].real))
    f.write('variable b221_im equal {:.9f} \n'.format(B1_inv[1,0].imag))
    f.write('variable b222_re equal {:.9f} \n'.format(B1_inv[1,1].real))
    f.write('variable b222_im equal {:.9f} \n'.format(B1_inv[1,1].imag))
    f.write('variable b223_re equal {:.9f} \n'.format(B1_inv[1,2].real))
    f.write('variable b223_im equal {:.9f} \n'.format(B1_inv[1,2].imag))
    f.write('variable b231_re equal {:.9f} \n'.format(B1_inv[2,0].real))
    f.write('variable b231_im equal {:.9f} \n'.format(B1_inv[2,0].imag))
    f.write('variable b232_re equal {:.9f} \n'.format(B1_inv[2,1].real))
    f.write('variable b232_im equal {:.9f} \n'.format(B1_inv[2,1].imag))
    f.write('variable b233_re equal {:.9f} \n'.format(B1_inv[2,2].real))
    f.write('variable b233_im equal {:.9f} \n'.format(B1_inv[2,2].imag))
    
    f.write('\n')
    f.write('# Interfacial parameters\n')
    
    # S_breve-matrix
    f.write('# S_breve-matrix\n')
    
    f.write('variable s_11 equal {:.9f} \n'.format(S_brev[0,0]))
    f.write('variable s_12 equal {:.9f} \n'.format(S_brev[0,1]))
    f.write('variable s_13 equal {:.9f} \n'.format(S_brev[0,2]))
    f.write('variable s_21 equal {:.9f} \n'.format(S_brev[1,0]))
    f.write('variable s_22 equal {:.9f} \n'.format(S_brev[1,1]))
    f.write('variable s_23 equal {:.9f} \n'.format(S_brev[1,2]))
    f.write('variable s_31 equal {:.9f} \n'.format(S_brev[2,0]))
    f.write('variable s_32 equal {:.9f} \n'.format(S_brev[2,1]))
    f.write('variable s_33 equal {:.9f} \n'.format(S_brev[2,2]))
    
    # Phi_inv-matrix
    f.write('# Phi_inv-matrix\n')
    f.write('variable Phi_11 equal {:.9f} \n'.format(Phi_inv[0,0]))
    f.write('variable Phi_12 equal {:.9f} \n'.format(Phi_inv[0,1]))
    f.write('variable Phi_13 equal {:.9f} \n'.format(Phi_inv[0,2]))
    f.write('variable Phi_21 equal {:.9f} \n'.format(Phi_inv[1,0]))
    f.write('variable Phi_22 equal {:.9f} \n'.format(Phi_inv[1,1]))
    f.write('variable Phi_23 equal {:.9f} \n'.format(Phi_inv[1,2]))
    f.write('variable Phi_31 equal {:.9f} \n'.format(Phi_inv[2,0]))
    f.write('variable Phi_32 equal {:.9f} \n'.format(Phi_inv[2,1]))
    f.write('variable Phi_33 equal {:.9f} \n'.format(Phi_inv[2,2]))
    
    f.write('\n')
    f.write('variable E_22 equal {:.9f} \n'.format(E[1,1]))

else:
    raise Exception("Single_crytal flag was incorrectly specified!")
f.close()
