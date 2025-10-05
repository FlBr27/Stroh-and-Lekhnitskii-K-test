# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:19:59 2024

@author: Brunn

This function library contains functions that were taken from: 
    
Cheng et al., (2018). Aimsgb: An algorithm and open-source python library to
generate periodic grain boundary structures
Computational Materials Science 155 (2018) 92–103
https://doi.org/10.1016/j.commatsci.2018.08.029
and
https://github.com/ksyang2013/aimsgb

and 

Hadian et al., (2018). GB code: A grain boundary generation code.
Journal of Open Source Software, 3(29), 900.
https://doi.org/10.21105/joss.00900
and
https://github.com/oekosheri/GB_code/blob/master/gb_code/csl_generator.py

These functiond will be also cited individually below.
"""

import os
import numpy as np 
from numpy import linalg as la
import math
import warnings
from cxroots import Rectangle

###############################################################################
def is_positive_definite(M):
    """
    checks if M is positive definite (i.e. all eigenvalues are positive)
    """
    
    definite_check = np.all(la.eigvals(M) > 0)
    
    return definite_check

###############################################################################
def build_iso_stiff(mu_L, lambda_L):
    """
    Parameters
    ----------
    mu_L, lambda_L : float; Lame constants

    Returns
    -------
    C : 6 x 6 array
        stiffness matrix in Voigt notation
    """
    
    C_11 = lambda_L + 2*mu_L
    C_12 = lambda_L
    C_44 = mu_L

    C = np.array([[C_11, C_12,   C_12,   0,      0,      0],
                  [C_12, C_11,   C_12,   0,      0,      0],
                  [C_12, C_12,   C_11,   0,      0,      0],
                  [0,    0,      0,      C_44,   0,      0],
                  [0,    0,      0,      0,      C_44,   0],
                  [0,    0,      0,      0,      0,      C_44]])

    return C

###############################################################################
def build_cubic_stiff(C_11, C_12, C_44):
    """
    Parameters
    ----------
    C_11 : float
        Elastic constant C_11 of a bcc or fcc crystal
    C_12 : float
        Elastic constant C_12 of a bcc or fcc crystal
    C_44 : float
        Elastic constant C_44 of a bcc or fcc crystal

    Returns
    -------
    C : 6 x 6 array
        stiffness matrix in Voigt notation
    """

    C = np.array([[C_11, C_12,   C_12,   0,      0,      0],
                  [C_12, C_11,   C_12,   0,      0,      0],
                  [C_12, C_12,   C_11,   0,      0,      0],
                  [0,    0,      0,      C_44,   0,      0],
                  [0,    0,      0,      0,      C_44,   0],
                  [0,    0,      0,      0,      0,      C_44]])

    return C

###############################################################################
def build_hex_stiff(C_11, C_12, C_13, C_33, C_44):
    """
    Parameters
    ----------
    C_11 : float
        Elastic constant C_11 of an hcp crystal
    C_12 : float
        Elastic constant C_12 of an hcp crystal
    C_13 : float
        Elastic constant C_13 of an hcp crystal
    C_44 : float
        Elastic constant C_44 of an hcp crystal

    Returns
    -------
    C : 6 x 6 array
        stiffness matrix in Voigt notation
    """

    C = np.array([[C_11, C_12,   C_13,   0,      0,      0],
                  [C_12, C_11,   C_13,   0,      0,      0],
                  [C_13, C_13,   C_33,   0,      0,      0],
                  [0,    0,      0,      C_44,   0,      0],
                  [0,    0,      0,      0,      C_44,   0],
                  [0,    0,      0,      0,      0,      0.5*(C_11-C_12)]])

    return C

###############################################################################
def rotate_stiff_mat(C,x,y,z):
    """
    Rotates a stiffness matrix in VOIGT NOTATION into a crystallographic
    orientation, specified by the vectors x, y, and z.

    Parameters
    ----------
    C : 6 x 6 array
        Elasticity matrix in crystallographic standard coordinate system
    x : 3 x 1 array
        desired x-direction
    y : 3 x 1 array
        desired y-direction
    z : 3 x 1 array
        desired z-direction

    Returns
    -------
    C_rot : 6 x 6 array
        Elasticity matrix in rotated coordinate system
    """

    # normalized vectors
    x_n = (1./la.norm(x))*x
    y_n = (1./la.norm(y))*y
    z_n = (1./la.norm(z))*z

    # sanity check 1
    eps = 1.0e-12
    sc1 = la.norm(np.cross(x_n, y_n) - z_n)
    sc2 = la.norm(np.cross(y_n, z_n) - x_n)

    if sc1 > eps or sc2 > eps:
        raise Exception("The specified set of vectors is not an orthogonal "
                        "right-handed coordinate system!")

    # standard orientations of the crystal
    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    # 3x3 rotation matrix
    R = np.array([[np.dot(x_n,e1), np.dot(x_n,e2), np.dot(x_n,e3)],
                  [np.dot(y_n,e1), np.dot(y_n,e2), np.dot(y_n,e3)],
                  [np.dot(z_n,e1), np.dot(z_n,e2), np.dot(z_n,e3)]])

    # sanity check 2
    if not -eps < la.det(R) - 1. < eps:
        raise Exception("The determinant of the 3x3 rotation matrix is not 1!")

    # 6x6 rotation matrix according to Aud, p.74

    R6_11 = np.array([[R[0,0]**2, R[0,1]**2, R[0,2]**2],
                      [R[1,0]**2, R[1,1]**2, R[1,2]**2],
                      [R[2,0]**2, R[2,1]**2, R[2,2]**2]])

    R6_12 = np.array([[R[0,1]*R[0,2], R[0,2]*R[0,0], R[0,0]*R[0,1]],
                      [R[1,1]*R[1,2], R[1,2]*R[1,0], R[1,0]*R[1,1]],
                      [R[2,1]*R[2,2], R[2,2]*R[2,0], R[2,0]*R[2,1]]])

    R6_21 = np.array([[R[1,0]*R[2,0], R[1,1]*R[2,1], R[1,2]*R[2,2]],
                      [R[2,0]*R[0,0], R[2,1]*R[0,1], R[2,2]*R[0,2]],
                      [R[0,0]*R[1,0], R[0,1]*R[1,1], R[0,2]*R[1,2]]])

    R6_22 = np.array([[R[1,1]*R[2,2]+R[1,2]*R[2,1], R[1,0]*R[2,2]+R[1,2]*R[2,0],
                      R[1,1]*R[2,0]+R[1,0]*R[2,1]],

                      [R[0,1]*R[2,2]+R[0,2]*R[2,1], R[0,2]*R[2,0]+R[0,0]*R[2,2],
                                        R[0,0]*R[2,1]+R[0,1]*R[2,0]],

                      [R[0,1]*R[1,2]+R[0,2]*R[1,1], R[0,2]*R[1,0]+R[0,0]*R[1,2],
                                        R[0,0]*R[1,1]+R[0,1]*R[1,0]]])

    R6 = np.empty((6,6))

    R6[0:3, 0:3] = R6_11
    R6[0:3, 3:6] = 2*R6_12
    R6[3:6, 0:3] = R6_21
    R6[3:6, 3:6] = R6_22

    # sanity check 3
    if not -eps < la.det(R6) - 1. < eps:
        raise Exception("The determinant of the 6x6 rotation matrix is not 1!")

    C_rot = R6 @ C @ R6.T

    return C_rot

###############################################################################
def ang(a, b):
    """
    returns the cos(angle) between two vectors. [Hadian et al.]
    """
    ang = np.round(np.dot(a, b)/la.norm(a)/la.norm(b), 7)
    return abs(ang)

###############################################################################
def is_integer(a, tol=1e-5):
    """
    Check whether the number is integer; [Cheng et al.]
    """
    return la.norm(np.abs(a - np.round(a))) < tol

###############################################################################
def SymmEquivalent(arr):
    """
    returns cubic symmetric eqivalents of the given 2 dimensional vector.
    [Hadian et al.]
    """
    Sym = np.zeros([24, 3, 3])
    Sym[0, :] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Sym[1, :] = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    Sym[2, :] = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
    Sym[3, :] = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    Sym[4, :] = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
    Sym[5, :] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    Sym[6, :] = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
    Sym[7, :] = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    Sym[8, :] = [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
    Sym[9, :] = [[-1, 0, 0], [0, 0, 1], [0, 1, 0]]
    Sym[10, :] = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    Sym[11, :] = [[1, 0, 0], [0, 0, 1], [0, -1, 0]]
    Sym[12, :] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    Sym[13, :] = [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]
    Sym[14, :] = [[0, -1, 0], [0, 0, 1], [-1, 0, 0]]
    Sym[15, :] = [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
    Sym[16, :] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    Sym[17, :] = [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
    Sym[18, :] = [[0, 0, -1], [1, 0, 0], [0, -1, 0]]
    Sym[19, :] = [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
    Sym[20, :] = [[0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    Sym[21, :] = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
    Sym[22, :] = [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    Sym[23, :] = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

    arr = np.atleast_2d(arr)
    Result = []
    for i in range(len(Sym)):
        for j in range(len(arr)):
            Result.append(np.dot(Sym[i, :], arr[j]))
    Result = np.array(Result)
    return np.unique(Result, axis=0)

###############################################################################
def stgb_struc_search(y1, z1, R):
    """
    Calculates from the grain boundary normal vector y1, the common
    rotation axis z1 and the rotation matrix R the remaining vectors x1,
    x2, y2, z2 to fully specify both grains and returns also a flag that
    indicates if the determined vectors span a symmetric tilt grain boundary.
    """

    eps = 1.0e-6

    if np.abs(np.dot(y1,z1)) > eps:
        raise Exception("Specified y- and z-vector are not orthogonal!")

    x1 = np.cross(y1, z1)

    # simplify x1 if possible
    gcd_x1 = math.gcd(math.gcd(x1[0], x1[1]),x1[2])
    x1 = (x1/gcd_x1).astype(int)

    ### calculate the reamaing quantities ###
    x2_re = np.matmul(R,x1)
    x2 = np.round(x2_re).astype(int)
    y2_re = np.matmul(R,y1)
    y2 = np.round(y2_re).astype(int)
    z2 = z1

    ### checks ###

    # check if the rotated vector entries are integers
    sc1 = is_integer(x2_re, tol=1e-4)
    sc2 = is_integer(y2_re, tol=1e-4)

    # check if the the mean plane is a mirror plane of the cubic system
    # The theory behind this is explained e.g. in:
    # An et al. 1982: phys. stat. sol. (b) 114, 349 (1982)
    # The implementation follows the one of Hadian et al.
    # see Journal of Open Source Software, 3(29), 900.

    y_mean = (y1 + y2)/2
    sc3 = False

    # Mirrorplanes cubic symmetry
    MP = np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1],
                   [1, 1, 0],
                   [0, 1, 1],
                   [1, 0, 1],
                   ], dtype='float')

    MP_symm = SymmEquivalent(MP)


    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        for j in range(len(MP_symm)):
            if 1 - ang(y_mean, MP_symm[j]) < eps:
                sc3 = True
                break
    # Runtime warnings are filtered here, since csl.ang = nan doesn't harm the 
    # calculations but produces a warning
    
    flag = sc1 and sc2 and sc3

    return flag, x1, x2, y2, z2

###############################################################################
def find_stgb_struc(phi_deg, y1, z1):

    """
    Takes a symmetrically tilt grain boundary as an argument and returns
    all lattice vectors of both grains.
    The grain boundaries are expected to be specified in the CSL convention:
    Sigma X (y1) phi° / [z1] (e.g. Sigma 3 (1 1 1) 70.5° [1 -1 0])

    Parameters
    ----------
    phi_deg : float
        Rotation angle in degree
    y1 : 3 x 1 array
        Grain boundary normal w.r.t. grain 1
    z1 : 3 x 1 array
        Rotation axis common to both grains


    Returns
    -------
    x1 : 3 x 1 array
        x-direction of grain 1
    x2 : 3 x 1 array
        x-direction of grain 2
    y2 : 3 x 1 array
        y-direction of grain 2
    z2 : 3 x 1 array
        z-direction of grain 2
    """

    eps = 1.0e-12

    phi = phi_deg*(math.pi/180)
    z1_n = (1/np.linalg.norm(z1))*z1

    ### Rotation matrix for rotating coordinate system 1 clockwise into 
    ### coordinate system 2.
    # In this way, all vectors that are expressed in coordinate system 1
    # get expressed in coordinate system 2 and the two coordinate systems
    # can be brought in coincidence upon counter-clockwise (positive) rotation
    # of coordinate system 2
    
    sin_th = math.sin(phi)
    cos_th = math.cos(phi)

    R = np.empty((3 ,3))

    R[0,0] = cos_th + z1_n[0]**2*(1. - cos_th)
    R[0,1] = z1_n[0]*z1_n[1]*(1. - cos_th) - z1_n[2]*sin_th
    R[0,2] = z1_n[0]*z1_n[2]*(1. - cos_th) + z1_n[1]*sin_th

    R[1,0] = z1_n[1]*z1_n[0]*(1. - cos_th) + z1_n[2]*sin_th
    R[1,1] = cos_th + z1_n[1]**2*(1. - cos_th)
    R[1,2] = z1_n[1]*z1_n[2]*(1. - cos_th) - z1_n[0]*sin_th

    R[2,0] = z1_n[2]*z1_n[0]*(1. - cos_th) - z1_n[1]*sin_th
    R[2,1] = z1_n[2]*z1_n[1]*(1. - cos_th) + z1_n[0]*sin_th
    R[2,2] = cos_th + z1_n[2]**2*(1. - cos_th)


    ### do some sanity checks ###

    # determinant of R
    if not np.abs(np.linalg.det(R)) - 1. < eps:
        raise Exception("Determinant of the rotation matrix not equal to 1!")

    # R*R^T
    if not np.abs(np.trace(np.matmul(R, R.T) - np.eye(3))) < eps:
        raise Exception("There is a problem with the rotation matrix!")

    # invariance of z1 under R
    if not np.abs(np.linalg.norm(np.matmul(R,z1) - z1)) < eps:
        raise Exception("The rotation axis is not invariant w.r.t."
                        "the rotation matrix!")

    ### calculate the reamaing quantities ###
    flag, x1, x2, y2, z2 = stgb_struc_search(y1, z1, R)

    if not flag:
        raise Exception("Unable to find the required vectors for the "
                        "specified STGB!")

    return x1, x2, y2, z2

###############################################################################
def Qu_Bassani_check(C_el1, C_el2, silent=False):
    """
    | Checks if a grain boundary will yield oscillatory behaviour.
    | For silent=False (default) a message is printed to the
      screen, telling the result of the check.

    Parameters
    ----------
    C_el1 : 6 x 6 array
        Stiffness matrix in Voigt notation of Grain 1
    C_el2 : 6 x 6 array
        Stiffness matrix in Voigt notation of Grain 2

    Yields
    ------
    QB_flag : 'True' if GB will show oscillatory behavior and 'False' otherwise
    """

    # See Qu & Bassani, J. Mech. Phys. Solids Vol. 37. No. 4. pp. 417-433. 1989
    # for the theory
    
    eps = 1.0e-12
    check_1 = False
    check_2 = False

    if C_el1[0,3] > eps or C_el1[0,4] > eps or C_el1[1,3] > eps \
        or C_el1[1,4] > eps or C_el1[2,3] > eps or C_el1[2,4] > eps \
            or C_el1[3,5] > eps or C_el1[4,5] > eps:
                check_1 = True

    if C_el2[0,3] > eps or C_el2[0,4] > eps or C_el2[1,3] > eps \
        or C_el2[1,4] > eps or C_el2[2,3] > eps or C_el2[2,4] > eps \
            or C_el2[3,5] > eps or C_el2[4,5] > eps:
                check_2 = True

    if (check_1 or check_2):
        QB_flag = True
        if not silent:
            print("The specified grain boundary will yield oscillatory "
                  "behaviour!")

    else:
        QB_flag = False
        if not silent:
            print("The specified grain boundary will NOT yield "
                  "oscillatory behaviour!")

    return QB_flag

###############################################################################
def Stroh6(C_el):
    """
    | Takes a (rotated) elasticity Matrix in Voigt notation and computes the
      eigenvalues p_i and the matrices A and B according to the 6th order
      Stroh formalism.
    | For details see Chapter 5 of  T.C.T. Ting, Anisotropic Elasticity -
      Theory and Applications, Oxford University Press, 1996

    Parameters
    ----------
    C_el : 6 x 6 float array; rotated stiffness matrix in Voigt notation

    """

    # all equations without further reference are from:
    # T.C.T. Ting, Anisotropic Elasticity - Theory and Applications,
    # Oxford University Press, 1996

    eps1 = 1.0e-12
    eps2 = 1.0e-6

    # Q, R and T matrices (Eq. 5.2-1)
    Q = np.array([[C_el[0,0], C_el[0,5], C_el[0,4]],
                  [C_el[0,5], C_el[5,5], C_el[4,5]],
                  [C_el[0,4], C_el[4,5], C_el[4,4]]])

    R = np.array([[C_el[0,5], C_el[0,1], C_el[0,3]],
                  [C_el[5,5], C_el[1,5], C_el[3,5]],
                  [C_el[4,5], C_el[1,4], C_el[3,4]]])

    T = np.array([[C_el[5,5], C_el[1,5], C_el[3,5]],
                  [C_el[1,5], C_el[1,1], C_el[1,3]],
                  [C_el[3,5], C_el[1,3], C_el[3,3]]])

    # fundamental elasticity matrix (Eq. 5.5-4 and 5.5-5)
    T_inv = la.inv(T)
    N1 = -T_inv @ R.T
    N2 = T_inv
    N3 = R @ T_inv @ R.T - Q

    N = np.zeros([6,6])

    N[0:3,0:3] = N1
    N[0:3,3:6] = N2
    N[3:6,0:3] = N3
    N[3:6,3:6] = N1.T

    eigVal, eigVec = la.eig(N)

    mask = np.imag(eigVal) > 0
    p = eigVal[mask]
    xi = eigVec[:,mask]

    # check if the system is degenerated
    dgc1 = np.abs(p[0] - p[1]) < eps2
    dgc2 = np.abs(p[1] - p[2]) < eps2
    dgc3 = np.abs(p[0] - p[2]) < eps2

    if dgc1 or dgc2 or dgc3:
       raise Exception("The specified elasticity matrix is degenerated "
                        "in terms of the Stroh formalism!")

    # normalization of xi
    I3 = np.eye(3)
    I6 = np.array([[0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 0, 1],
                   [1, 0, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0]])

    eta = I6 @ xi

    normf = np.sqrt(eta.T @ xi)

    xi_n = np.array([xi[:,0]/normf[0,0], xi[:,1]/normf[1,1],\
                     xi[:,2]/normf[2,2]]).T
    eta_n = I6 @ xi_n

    # verify normalization
    normcheck = eta_n.T @ xi_n
    nc1 = np.real(normcheck[0,0]) - 1.0 < eps1
    nc2 = np.real(normcheck[1,1]) - 1.0 < eps1
    nc3 = np.real(normcheck[2,2]) - 1.0 < eps1
    nc4 = np.imag(normcheck[0,0]) < eps1
    nc5 = np.imag(normcheck[1,1]) < eps1
    nc6 = np.imag(normcheck[2,2]) < eps1
    nc7 = np.abs(la.norm(normcheck) - la.norm(I3)) < eps1

    if not (nc1 and nc2 and nc3 and nc4 and nc5 and nc6 and nc7):
        raise Exception("Stroh eigenvectors are not properly normalized!")

    # A and B matrices (Eq. 5.5-4)
    A = xi_n[0:3,:]
    B = xi_n[3:6,:]

    return p, A, B

###############################################################################
def Stroh6_Coeff(p1, A1, B1, p2, A2, B2, a):
    """
    | Calculates the required parameters to determine the displacement and
      stress fields according to F. Brunner's solutions to the 6th
      order Stroh problem for interfacial cracks.
    | In the following region 1 corresponds to the upper half-body (y >= 0)
      and region 2 corresponds to the lower half-body (y < 0).

    Parameters
    ----------
    | p1, p2: Stroh eigenvalues of the corresponding region [-];
      3x1 float arrays;
    | A1, B1, A2, B2: Stroh matrices of the corresponding region;
      A_i [1/GPa], B_i [-]; 3x3 float arrays;
    | a: crack length [A]; float;

    Returns
    -------
    | Param: dictonary that contains all relevant parameters for the
      calculation of the field quantities that are not depending on the
      coordinates.
    """

    # all equations without further reference are from:
    # T.C.T. Ting, Anisotropic Elasticity - Theory and Applications,
    # Oxford University Press, 1996

    eps = 1.0e-12

    # inverse impedance tensors (Eq. 5.6-1)
    M_inv1 = 1j*A1 @ la.inv(B1)
    M_inv2 = 1j*A2 @ la.inv(B2)

    # bimaterial matrix
    # cf. C. Hwu, Anisotropic Elastic Plates, Springer, 2010, Eq. 7.11b
    M_st = M_inv1 + np.conjugate(M_inv2)

    D = np.real(M_st)
    W = -np.imag(M_st)

    # check if D is symmetric and W anti-symmetric as required
    symcheck1 = la.norm(D - D.T) < eps
    symcheck2 = la.norm(W + W.T) < eps

    if not (symcheck1 and symcheck2):
        raise Exception("An Error occured in the calculation of the "
                        "D and W matrices!")

    # S_breve, beta and gamma (Eq. 11.6-2 - 11.6-4)
    pi = np.pi
    I = np.eye(3)
    S_brev = la.inv(D) @ W

    # Distinction between oscillatory and non-oscillatory case
    beta_tilde = np.trace(S_brev @ S_brev)

    if np.abs(beta_tilde) <= eps:

        beta = 0.0
        gamma = 0.0
        S_brev = np.zeros([3,3])
        Phi_inv = I
        oscill_flag = False

    else:

        beta = (-0.5*beta_tilde)**(1/2)

        # validate beta (Eq. 11.6-3)
        beta_chk1 = np.isreal(beta)
        beta_chk2 = 0.0 <= beta < 1.0

        if not (beta_chk1 and beta_chk2):
            raise Exception("Parameter beta was calculated incorrectly!")

        gamma = 1/pi*np.arctanh(beta)

        S_brev3 = S_brev @ S_brev @ S_brev
        Phi = I + 2*np.arctanh(beta)/(pi*beta**3)*S_brev3
        Phi_inv = la.inv(Phi)
        oscill_flag = True


    p1_mat = np.array([[p1[0], 0, 0],
                       [0, p1[1], 0],
                       [0, 0, p1[2]]])

    p2_mat = np.array([[p2[0], 0, 0],
                       [0, p2[1], 0],
                       [0, 0, p2[2]]])

    Param = {
        'p1_mat': p1_mat,
        'A1': A1,
        'B1': B1,
        'p2_mat': p2_mat,
        'A2': A2,
        'B2': B2,
        'S_brev': S_brev,
        'gamma': gamma,
        'beta': beta,
        'Phi_inv': Phi_inv,
        'a': a,
        'oscill_flag': oscill_flag
        }

    return Param

###############################################################################
def LES(C_el, mode='lekh4'):
    """
    Lekhnitskii-Eshelby-Stroh formalism to obtain the matrices A and B,
    and the eigenvalues mu_j (p_j) that are used in the Stroh formalism,
    or the inverse impedance matrix M^-1, the reduced compliance Matrix S_red
    and mu_j, that are used in the Lekhnitskii formalism in terms of the
    4th-order Lekhnitskii formalism.

    Parameters
    ----------
    | C_el : 6 x 6 float array; (rotated) stiffness matrix in Voigt notation
    | mode :
    | 'stroh4': A, B, mu of the 4th order plane strain formalism are returned
    | 'lekh4': M^-1, S_red, mu of the 4th order plane strain formalism are
      returned

    Returns
    -------
    | A: 3 x 3 float array
    | B: 3 x 3 float array
    | mu: 3 x 1 float array
    | OR
    | M_inv: 3 x 3 float array
    | S_red: 6 x 6 float array
    | mu: 3 x 1 float array
    """

    modes = ['stroh4', 'lekh4']

    if mode not in modes:
        raise Exception("Invalid mode specification!")
        
    eps = 1.0e-12
        
    # Check if the (classical) plain strain assumption holds
    if not (np.abs(C_el[0,3]) < eps and np.abs(C_el[0,4]) < eps and
            np.abs(C_el[1,3]) < eps and np.abs(C_el[1,4]) < eps and
            np.abs(C_el[2,3]) < eps and np.abs(C_el[2,4]) < eps and
            np.abs(C_el[3,5]) < eps and np.abs(C_el[4,5]) < eps):
    
                raise Exception("The specified elasticity matrix violates "
                                "the classical plane strain assumption!")

    # compute compliance matrix
    S_el =  la.inv(C_el)
    S_red = np.empty([6,6])

    # compute reduced compliance matrix
    for i1 in range(6):
        for j1 in range(6):
            S_red[i1,j1] = S_el[i1,j1] - (S_el[i1,2]*S_el[2,j1])/S_el[2,2]


    # finding the roots of the characteristic equations 
    # cf. Ting - Anisotropic Elasticity, p. 125
   
    # setting up the 4th-order polynomial by it's coefficients
    # [p^4, p^3, p^2, p^1, p^0]
    coeff = [S_red[0,0], -2*S_red[0,5], (2*S_red[0,1] + S_red[5,5]), 
            -2*S_red[1,5], S_red[1,1]]

    # finding the first two roots
    mu = (np.roots(coeff)).tolist()
    
    # find third root
    mu3 = (S_red[3,4] + 1j*np.sqrt(S_red[3,3]*S_red[4,4] -
                                  S_red[3,4]**2))/S_red[4,4]

    ### Check roots ###
    # NOTE: The treatment of repeated roots used here follows the approximate
    # approach as suggested by C. Hwu in "Anisotropic Elastic Plates", 2010,
    # p. 86. The proper treatment as suggested in Section 3.5 would be better,
    # however, it is also more involved.
    if len(mu) == 4:

        # keep only roots with Im(mu) > 0
        mu = [x for x in mu if np.imag(x) >= 0.0]

    elif len(mu) == 2:
        # if only one pair of roots was found, both roots are used but the one
        # with Im(mu) < 0 is conjugated. In this way, two slightly different
        # roots are obtained.

        print("Warning: the material is degenerated w.r.t. the "
              "Lekhnitskii formalism!")

        if np.imag(mu[0]) < 0:
            mu[0] = np.conj(mu[0])

        if np.imag(mu[1]) < 0:
            mu[1] = np.conj(mu[1])

    else: raise Exception("An error occurred during the root determination!")

    mu.append(mu3)
    
    mu = np.array(mu,dtype=complex)
    
    ### calculate the A and B matrix
    # see Suo, Proc. R. Soc. Lond. A 427, 331-358, 1990 (A = A and B = L there)

    B = np.array([[-mu[0], -mu[1], 0],
                  [  1,      1,    0],
                  [  0,      0,   -1]])

    A = np.empty([3,3],dtype=complex)
    
    A[0,0] = S_red[0,0]*mu[0]**2 + S_red[0,1] - S_red[0,5]*mu[0]

    A[0,1] = S_red[0,0]*mu[1]**2 + S_red[0,1] - S_red[0,5]*mu[1] 

    A[0,2] = 0.0

    A[1,0] = S_red[1,0]*mu[0] + S_red[1,1]/mu[0] - S_red[1,5] 

    A[1,1] = S_red[1,0]*mu[1] + S_red[1,1]/mu[1] - S_red[1,5] 

    A[1,2] = 0.0

    A[2,0] = 0.0

    A[2,1] = 0.0

    A[2,2] = S_red[3,4] - S_red[3,3]/mu[2]

    # calculate the impedance matrix M^-1 from A and B
    # see Wu & Curtin, Acta Materialia 88, 1-12, 2015
    M_inv = 1j*A @ la.inv(B)

    if  mode == 'stroh4':
        return np.round(A,15), np.round(B,15), np.round(mu,15)
    elif mode == 'lekh4':
        return np.round(M_inv,15), np.round(S_red,15), np.round(mu,15)
    # rounds to get rid of values < O(10^-15)

###############################################################################
def GB_InpGen_Lekh_4ps(C_el1, C_el2, write='None',
                       outpf_name='GB_InpGen_Lekh_4ps_out',ex_hand='on'):
    """
    Takes the two rotated stiffness matrices C_el1 and C_el2 of a grain
    boundary in Voigt notation and generates the input parameters for LAMMPS
    K-tests of non-oscillatory grain boundaries according to the 4th order
    plane strain Lekhnitskii formalism.

    Parameters
    ----------
    | C_el1, C_el2:  6 x 6 float array; rotated stiffness matrices in Voigt
      notation
    |
    | write: string;
    | 'None': No additional output file is created; default option.
    | 'Info': The parameters are saved internally to a .txt file in
       the 'output_files' folder.
    | 'LAMMPS': The parameters are saved to a .in file in the
      'K-Test/include_inputs/Input_Parameters_no' folder such that it can
      be directly written by LAMMPS scripts
    |
    | outpf_name: string; Name of the output file if created. The default is
      'GB_InpGen_Lekh_4ps_out'.
    | ex_hand: 'on' or 'off'; determines if exceptions are raised if e.g. the
      bimaterial matrix is not real. Should only be turned off when you know
      what you are doing! The default is 'on'.

    Returns
    -------
    | coeff : dictonary; Contains all relevant parameters.
    """
    writes = ['None', 'Info', 'LAMMPS']

    if write not in writes:
        raise Exception("Invalid write style!")

    eps = 1.0e-12

    # get the two impedance matrices
    M_inv1, S_red1, mu1 = LES(C_el1, mode='lekh4')
    M_inv2, S_red2, mu2 = LES(C_el2, mode='lekh4')

    # determine bimaterial matrix
    M_st = M_inv1 + np.conjugate(M_inv2)

    # check if M_st is real (non-oscillatory)
    CM = np.sum(np.imag(M_st))
    if (not -eps < CM < eps) and (ex_hand == 'on'):
        raise Exception("The bimaterial matrix is not real!")

    D = np.real(M_st)

    # calculate the required parameters for LAMMPS
    p11 = S_red1[0,0]*mu1[0]**2 + S_red1[0,1] - S_red1[0,5]*mu1[0]
    p12 = S_red1[0,0]*mu1[1]**2 + S_red1[0,1] - S_red1[0,5]*mu1[1]
    q11 = (S_red1[0,1]*mu1[0]**2 + S_red1[1,1] - S_red1[1,5]*mu1[0])/mu1[0]
    q12 = (S_red1[0,1]*mu1[1]**2 + S_red1[1,1] - S_red1[1,5]*mu1[1])/mu1[1]

    p21 = S_red2[0,0]*mu2[0]**2 + S_red2[0,1] - S_red2[0,5]*mu2[0]
    p22 = S_red2[0,0]*mu2[1]**2 + S_red2[0,1] - S_red2[0,5]*mu2[1]
    q21 = (S_red2[0,1]*mu2[0]**2 + S_red2[1,1] - S_red2[1,5]*mu2[0])/mu2[0]
    q22 = (S_red2[0,1]*mu2[1]**2 + S_red2[1,1] - S_red2[1,5]*mu2[1])/mu2[1]

    coeff = {'s11': mu1[0],
             's12': mu1[1],
             'p11': p11,
             'p12': p12,
             'q11': q11,
             'q12': q12,
             's21': mu2[0],
             's22': mu2[1],
             'p21': p21,
             'p22': p22,
             'q21': q21,
             'q22': q22,
             'D22': D[1,1]}

    if write == 'Info':

        # writing parameters to file       
        outpf_name = './outputs/' + outpf_name + ".txt"
        os.makedirs(os.path.dirname(outpf_name), exist_ok=True)
        f = open(outpf_name, 'w')

        f.write('s11 = {:.7f}, {:.7f}i \n'.format(mu1[0].real, mu1[0].imag))
        f.write('s12 = {:.7f}, {:.7f}i \n'.format(mu1[1].real, mu1[1].imag))
        f.write('p11 = {:.7f}, {:.7f}i \n'.format(p11.real, p11.imag))
        f.write('p12 = {:.7f}, {:.7f}i \n'.format(p12.real, p12.imag))
        f.write('q11 = {:.7f}, {:.7f}i \n'.format(q11.real, q11.imag))
        f.write('q12 = {:.7f}, {:.7f}i \n'.format(q12.real, q12.imag))
        f.write('\n')
        f.write('s21 = {:.7f}, {:.7f}i \n'.format(mu2[0].real, mu2[0].imag))
        f.write('s22 = {:.7f}, {:.7f}i \n'.format(mu2[1].real, mu2[1].imag))
        f.write('p21 = {:.7f}, {:.7f}i \n'.format(p21.real, p21.imag))
        f.write('p22 = {:.7f}, {:.7f}i \n'.format(p22.real, p22.imag))
        f.write('q21 = {:.7f}, {:.7f}i \n'.format(q21.real, q21.imag))
        f.write('q22 = {:.7f}, {:.7f}i \n'.format(q22.real, q22.imag))
        f.write('\n')
        f.write('D22 = {:.7f}'.format(D[1,1]))

        f.close()

    if write == 'LAMMPS':

        # writing parameters to file       
        outpf_name = './outputs/' + outpf_name
        os.makedirs(os.path.dirname(outpf_name), exist_ok=True)
        f = open(outpf_name, 'w')

        f.write('# Input parameters for the non-oscillatory 4th order ')
        f.write('Lekhnitskii K-test \n')
        f.write('\n')
        f.write('variable s11_real equal {:.9f} \n'.format(mu1[0].real))
        f.write('variable s11_imag equal {:.9f} \n'.format(mu1[0].imag))
        f.write('variable s12_real equal {:.9f} \n'.format(mu1[1].real))
        f.write('variable s12_imag equal {:.9f} \n'.format(mu1[1].imag))
        f.write('variable p11_real equal {:.9f} \n'.format(p11.real))
        f.write('variable p11_imag equal {:.9f} \n'.format(p11.imag))
        f.write('variable p12_real equal {:.9f} \n'.format(p12.real))
        f.write('variable p12_imag equal {:.9f} \n'.format(p12.imag))
        f.write('variable q11_real equal {:.9f} \n'.format(q11.real))
        f.write('variable q11_imag equal {:.9f} \n'.format(q11.imag))
        f.write('variable q12_real equal {:.9f} \n'.format(q12.real))
        f.write('variable q12_imag equal {:.9f} \n'.format(q12.imag))
        f.write('\n')
        f.write('variable s21_real equal {:.9f} \n'.format(mu2[0].real))
        f.write('variable s21_imag equal {:.9f} \n'.format(mu2[0].imag))
        f.write('variable s22_real equal {:.9f} \n'.format(mu2[1].real))
        f.write('variable s22_imag equal {:.9f} \n'.format(mu2[1].imag))
        f.write('variable p21_real equal {:.9f} \n'.format(p21.real))
        f.write('variable p21_imag equal {:.9f} \n'.format(p21.imag))
        f.write('variable p22_real equal {:.9f} \n'.format(p22.real))
        f.write('variable p22_imag equal {:.9f} \n'.format(p22.imag))
        f.write('variable q21_real equal {:.9f} \n'.format(q21.real))
        f.write('variable q21_imag equal {:.9f} \n'.format(q21.imag))
        f.write('variable q22_real equal {:.9f} \n'.format(q22.real))
        f.write('variable q22_imag equal {:.9f} \n'.format(q22.imag))
        f.write('\n')
        f.write('variable D22 equal {:.9f}'.format(D[1,1]))

        f.close()

    return coeff

###############################################################################
def slip_sys_bcc():

# The below array contains the possible slip systems of a bcc crystal
# the first array defines the slip plane normal
# the second array defines the slip direction

    slip_sys_bcc = np.array(
        # {110}
        [[[1, 1, 0], [1, -1, 1]],
         [[1, 1, 0], [-1, 1, 1]],
         [[1, -1, 0], [1, 1, -1]],
         [[1, -1, 0], [1, 1, 1]],
         [[1, 0, -1], [1, 1, 1]],
         [[1, 0, -1], [1, -1, 1]],
         [[1, 0, 1], [1, 1, -1]],
         [[1, 0, 1], [-1, 1, 1]],
         [[0, 1, 1], [1, 1, -1]],
         [[0, 1, 1], [1, -1, 1]],
         [[0, 1, -1], [1, 1, 1]],
         [[0, 1, -1], [-1, 1, 1]],
         # {112}
         [[1, 1, 2], [1, 1, -1]],
         [[1, 2, 1], [1, -1, 1]],
         [[2, 1, 1], [-1, 1, 1]],
         [[1, 1, -2], [1, 1, 1]],
         [[1, -2, 1], [1, 1, 1]],
         [[-2, 1, 1], [1, 1, 1]],
         [[1, -1, 2], [-1, 1, 1]],
         [[1, 2, -1], [-1, 1, 1]],
         [[2, 1, -1], [1, -1, 1]],
         [[-1, 1, 2], [1, -1, 1]],
         [[-1, 2, 1], [1, 1, -1]],
         [[2, -1, 1], [1, 1, -1]],
         # {123}
         [[1, 2, 3], [1, 1, -1]],
         [[1, 3, 2], [1, -1, 1]],
         [[3, 1, 2], [-1, 1, 1]],
         [[2, 1, 3], [1, 1, -1]],
         [[2, 3, 1], [1, -1, 1]],
         [[3, 2, 1], [-1, 1, 1]],
         [[-1, 2, 3], [1, -1, 1]],
         [[-1, 3, 2], [1, 1, -1]],
         [[3, -1, 2], [1, 1, -1]],
         [[2, -1, 3], [-1, 1, 1]],
         [[2, 3, -1], [-1, 1, 1]],
         [[3, 2, -1], [1, -1, 1]],
         [[1, -2, 3], [-1, 1, 1]],
         [[1, 3, -2], [-1, 1, 1]],
         [[3, 1, -2], [1, -1, 1]],
         [[-2, 1, 3], [1, -1, 1]],
         [[-2, 3, 1], [1, 1, -1]],
         [[3, -2, 1], [1, 1, -1]],
         [[1, 2, -3], [1, 1, 1]],
         [[1, -3, 2], [1, 1, 1]],
         [[-3, 1, 2], [1, 1, 1]],
         [[2, 1, -3], [1, 1, 1]],
         [[2, -3, 1], [1, 1, 1]],
         [[-3, 2, 1], [1, 1, 1]]])
    
    return slip_sys_bcc

###############################################################################