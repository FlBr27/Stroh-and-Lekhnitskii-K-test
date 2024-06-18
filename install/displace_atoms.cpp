// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*
log: 
Sep 21, 2022: add style {crackaniso} -> apply displacement BCs. solved from anisotropic fracture for crack in a single crystal.
Oct 23, 2023: add style {gbcrack} -> apply displacement BCs to model cracks along grain boundaries.
*/
#include "displace_atoms.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_park.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <complex>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::MY_2PI;

//*************************************************************************
// Start of modification for gbcrack.
// *** Required complex linear algebra tools ***

// complex matrix struct
struct CMAT {
    std::complex<double> cmat[3][3];
};

// complex vector struct
struct CVEC {
    std::complex<double> cvec[3];
};

// complex 3x3 matrix - 3x3 matrix multiplication
CMAT cMxM(CMAT A, CMAT B)
{
    CMAT C;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C.cmat[i][j] = A.cmat[i][0]*B.cmat[0][j] + A.cmat[i][1]*B.cmat[1][j] + A.cmat[i][2]*B.cmat[2][j];
        }
    }

    return C;
}

// complex 3x3 matrix - 3x1 vector multiplication
CVEC cMxV(CMAT A, CVEC v)
{
    CVEC u;

    for (int i = 0; i < 3; ++i) {
        u.cvec[i] = A.cmat[i][0]*v.cvec[0] + A.cmat[i][1]*v.cvec[1] + A.cmat[i][2]*v.cvec[2];
    }

    return u;
}

// element-wise multiplication of a complex scalar with a complex 3x3 matrix
CMAT csxM(std::complex<double> s, CMAT A)
{
    CMAT C;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C.cmat[i][j] = s*A.cmat[i][j];
        }
    }

    return C;
}

// element-wise addition of two complex 3x3 matrices
CMAT cMpM(CMAT A, CMAT B)
{
    CMAT C;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C.cmat[i][j] = A.cmat[i][j] + B.cmat[i][j];
        }
    }

    return C;
}

// Modification ends.
//*************************************************************************

enum{MOVE,RAMP,RANDOM,ROTATE,CRACKANISO,GBCRACK};

/* ---------------------------------------------------------------------- */

DisplaceAtoms::DisplaceAtoms(LAMMPS *_lmp) : Command(_lmp)
{
  mvec = nullptr;
}

/* ---------------------------------------------------------------------- */

DisplaceAtoms::~DisplaceAtoms()
{
  memory->destroy(mvec);
}

/* ---------------------------------------------------------------------- */

void DisplaceAtoms::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR,"Displace_atoms command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace_atoms command");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot displace_atoms after "
               "reading restart file with per-atom info");

  if (comm->me == 0) utils::logmesg(lmp,"Displacing atoms ...\n");

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_atoms group ID");
  groupbit = group->bitmask[igroup];

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Attempting to displace atoms in rigid bodies");

  int style = -1;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"random") == 0) style = RANDOM;
  else if (strcmp(arg[1],"rotate") == 0) style = ROTATE;
  //*******************************
  // This is where has been modified.
  // add a new style to displace crack based on anisotropic solution: crackaniso.
  else if (strcmp(arg[1],"crackaniso") == 0) style = CRACKANISO;
  else if (strcmp(arg[1],"gbcrack") == 0) style = GBCRACK;
  // Modify ends.
  //*******************************
  else error->all(FLERR,"Illegal displace_atoms command");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == RANDOM) options(narg-6,&arg[6]);
  else if (style == ROTATE) options(narg-9,&arg[9]);
  //*******************************
  // This is where has been modified for crack aniso.

  // Read inputs: K ${s1_real} ${s1_imag} ${s2_real} ${s2_imag} ${p1_real} ${p1_imag} 
  // ${p2_real} ${p2_imag} ${q1_real} ${q1_imag} ${q2_real} ${q2_imag} ${xtip} ${ytip}
  else if (style == CRACKANISO) options(narg-17,&arg[17]);
  //*******************************
  // Start of modification for gbcrack.
  
  // Read inputs (67 parameters in total): 
  // ${K_I} ${K_II} ${K_III} ${a} ${reg} ${xtip} ${ytip}
  
  // ${p1_re} ${p1_im} ${p2_re} ${p2_im} ${p3_re} ${p3_im} 

  // ${a11_re} ${a11_im} ${a12_re} ${a12_im} ${a13_re} ${a13_im}
  // ${a21_re} ${a21_im} ${a22_re} ${a22_im} ${a23_re} ${a23_im}
  // ${a31_re} ${a31_im} ${a32_re} ${a32_im} ${a33_re} ${a33_im}

  // ${b11_re} ${b11_im} ${b12_re} ${b12_im} ${b13_re} ${b13_im}
  // ${b21_re} ${b21_im} ${b22_re} ${b22_im} ${b23_re} ${b23_im}
  // ${b31_re} ${b31_im} ${b32_re} ${b32_im} ${b33_re} ${b33_im}

  // ${s11} ${s12} ${s13}
  // ${s21} ${s22} ${s23}
  // ${s31} ${s32} ${s33}

  // ${Phi11} ${Phi12} ${Phi13}
  // ${Phi21} ${Phi22} ${Phi23}
  // ${Phi31} ${Phi32} ${Phi33}
 
  else if (style == GBCRACK) options(narg-69,&arg[69]);
  // Modification ends.

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms by 3-vector or specified variable(s)

  if (style == MOVE) {
    move(0,arg[2],xscale);
    move(1,arg[3],yscale);
    move(2,arg[4],zscale);
  }

  // move atoms in ramped fashion

  if (style == RAMP) {

    int d_dim = 0;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*utils::numeric(FLERR,arg[3],false,lmp);
      d_hi = xscale*utils::numeric(FLERR,arg[4],false,lmp);
    } else if (d_dim == 1) {
      d_lo = yscale*utils::numeric(FLERR,arg[3],false,lmp);
      d_hi = yscale*utils::numeric(FLERR,arg[4],false,lmp);
    } else if (d_dim == 2) {
      d_lo = zscale*utils::numeric(FLERR,arg[3],false,lmp);
      d_hi = zscale*utils::numeric(FLERR,arg[4],false,lmp);
    }

    int coord_dim = 0;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*utils::numeric(FLERR,arg[6],false,lmp);
      coord_hi = xscale*utils::numeric(FLERR,arg[7],false,lmp);
    } else if (coord_dim == 1) {
      coord_lo = yscale*utils::numeric(FLERR,arg[6],false,lmp);
      coord_hi = yscale*utils::numeric(FLERR,arg[7],false,lmp);
    } else if (coord_dim == 2) {
      coord_lo = zscale*utils::numeric(FLERR,arg[6],false,lmp);
      coord_hi = zscale*utils::numeric(FLERR,arg[7],false,lmp);
    }

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double fraction,dramp;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
        fraction = MAX(fraction,0.0);
        fraction = MIN(fraction,1.0);
        dramp = d_lo + fraction*(d_hi - d_lo);
        x[i][d_dim] += dramp;
      }
    }
  }

  // move atoms randomly
  // makes atom result independent of what proc owns it via random->reset()

  if (style == RANDOM) {
    auto random = new RanPark(lmp,1);

    double dx = xscale*utils::numeric(FLERR,arg[2],false,lmp);
    double dy = yscale*utils::numeric(FLERR,arg[3],false,lmp);
    double dz = zscale*utils::numeric(FLERR,arg[4],false,lmp);
    int seed = utils::inumeric(FLERR,arg[5],false,lmp);
    if (seed <= 0) error->all(FLERR,"Illegal displace_atoms random command");

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        random->reset(seed,x[i]);
        x[i][0] += dx * 2.0*(random->uniform()-0.5);
        x[i][1] += dy * 2.0*(random->uniform()-0.5);
        x[i][2] += dz * 2.0*(random->uniform()-0.5);
      }
    }

    delete random;
  }

  // rotate atoms by right-hand rule by theta around R
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // R0 = runit = unit vector for R
  // D = X - P = vector from P to X
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(theta) + B sin(theta)

  if (style == ROTATE) {
    double theta_new;
    double axis[3],point[3],qrotate[4],qnew[4];
    double a[3],b[3],c[3],d[3],disp[3],runit[3];
    double *quat;

    int dim = domain->dimension;
    point[0] = xscale*utils::numeric(FLERR,arg[2],false,lmp);
    point[1] = yscale*utils::numeric(FLERR,arg[3],false,lmp);
    point[2] = zscale*utils::numeric(FLERR,arg[4],false,lmp);
    axis[0] = utils::numeric(FLERR,arg[5],false,lmp);
    axis[1] = utils::numeric(FLERR,arg[6],false,lmp);
    axis[2] = utils::numeric(FLERR,arg[7],false,lmp);
    double theta = utils::numeric(FLERR,arg[8],false,lmp);
    if (dim == 2 && (axis[0] != 0.0 || axis[1] != 0.0))
      error->all(FLERR,"Invalid displace_atoms rotate axis for 2d");

    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Zero length rotation vector with displace_atoms");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;

    double angle = DEG2RAD*theta;
    double cosine = cos(angle);
    double sine = sin(angle);

    double qcosine = cos(0.5*angle);
    double qsine = sin(0.5*angle);
    qrotate[0] = qcosine;
    qrotate[1] = runit[0]*qsine;
    qrotate[2] = runit[1]*qsine;
    qrotate[3] = runit[2]*qsine;

    double ddotr;

    // flags for additional orientation info stored by some atom styles

    int ellipsoid_flag = atom->ellipsoid_flag;
    int line_flag = atom->line_flag;
    int tri_flag = atom->tri_flag;
    int body_flag = atom->body_flag;

    int theta_flag = 0;
    int quat_flag = 0;
    if (line_flag) theta_flag = 1;
    if (ellipsoid_flag || tri_flag || body_flag) quat_flag = 1;

    // AtomVec pointers to retrieve per-atom storage of extra quantities

    auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>( atom->style_match("ellipsoid"));
    auto avec_line = dynamic_cast<AtomVecLine *>( atom->style_match("line"));
    auto avec_tri = dynamic_cast<AtomVecTri *>( atom->style_match("tri"));
    auto avec_body = dynamic_cast<AtomVecBody *>( atom->style_match("body"));

    double **x = atom->x;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;
    int *body = atom->body;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    imageint *image = atom->image;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        // unwrap coordinate and reset image flags accordingly
        domain->unmap(x[i],image[i]);
        image[i] = ((imageint) IMGMAX << IMG2BITS) |
          ((imageint) IMGMAX << IMGBITS) | IMGMAX;

        d[0] = x[i][0] - point[0];
        d[1] = x[i][1] - point[1];
        d[2] = x[i][2] - point[2];
        ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
        c[0] = ddotr*runit[0];
        c[1] = ddotr*runit[1];
        c[2] = ddotr*runit[2];
        a[0] = d[0] - c[0];
        a[1] = d[1] - c[1];
        a[2] = d[2] - c[2];
        b[0] = runit[1]*a[2] - runit[2]*a[1];
        b[1] = runit[2]*a[0] - runit[0]*a[2];
        b[2] = runit[0]*a[1] - runit[1]*a[0];
        disp[0] = a[0]*cosine  + b[0]*sine;
        disp[1] = a[1]*cosine  + b[1]*sine;
        disp[2] = a[2]*cosine  + b[2]*sine;
        x[i][0] = point[0] + c[0] + disp[0];
        x[i][1] = point[1] + c[1] + disp[1];
        if (dim == 3) x[i][2] = point[2] + c[2] + disp[2];

        // theta for lines

        if (theta_flag && line[i] >= 0.0) {
          theta_new = fmod(avec_line->bonus[line[i]].theta+angle,MY_2PI);
          avec_line->bonus[atom->line[i]].theta = theta_new;
        }

        // quats for ellipsoids, tris, and bodies

        if (quat_flag) {
          quat = nullptr;
          if (ellipsoid_flag && ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[ellipsoid[i]].quat;
          else if (tri_flag && tri[i] >= 0)
            quat = avec_tri->bonus[tri[i]].quat;
          else if (body_flag && body[i] >= 0)
            quat = avec_body->bonus[body[i]].quat;
          if (quat) {
            MathExtra::quatquat(qrotate,quat,qnew);
            quat[0] = qnew[0];
            quat[1] = qnew[1];
            quat[2] = qnew[2];
            quat[3] = qnew[3];
          }
        }
      }
    }
  }

  //*******************************
  // This is where has been modified. 
  // Add a new style of displacement: crackaniso

  if (style == CRACKANISO) {
    double K1 = utils::numeric(FLERR,arg[2],false,lmp);
    double Ka = utils::numeric(FLERR,arg[3],false,lmp);
    double Kb = utils::numeric(FLERR,arg[4],false,lmp);
    double Kc = utils::numeric(FLERR,arg[5],false,lmp);
    double Kd = utils::numeric(FLERR,arg[6],false,lmp);
    double Ke = utils::numeric(FLERR,arg[7],false,lmp);
    double Kf = utils::numeric(FLERR,arg[8],false,lmp);
    double Kg = utils::numeric(FLERR,arg[9],false,lmp);
    double Kh = utils::numeric(FLERR,arg[10],false,lmp);
    double Kj = utils::numeric(FLERR,arg[11],false,lmp);
    double Kk = utils::numeric(FLERR,arg[12],false,lmp);
    double Kl = utils::numeric(FLERR,arg[13],false,lmp);
    double Km = utils::numeric(FLERR,arg[14],false,lmp);
    double KXc = utils::numeric(FLERR,arg[15],false,lmp);
    double KYc = utils::numeric(FLERR,arg[16],false,lmp);
    double Kpi = 3.1415926535898;

    // Each processor loop over the nlocal atoms it owns:
      double **x = atom->x;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) { // if atoms belong to the group specified in the command


      // calc delx, dely, delz        <------------------
        double x1 = x[i][0]-KXc;
        double x2 = x[i][1]-KYc;
        double r = sqrt(x1*x1+x2*x2);
        double teta = atan2(x2,x1);

        double sinteta = sin(teta);
        double costeta = cos(teta);
        double coef = K1*sqrt(2*r/Kpi);      

        std::complex<double> sinteta1 = sinteta;
        std::complex<double> costeta1 = costeta;
  
        std::complex<double> s1(Ka,Kb);
        std::complex<double> s2(Kc,Kd);
        std::complex<double> p1(Ke,Kf);
        std::complex<double> p2(Kg,Kh);
        std::complex<double> q1(Kj,Kk);
        std::complex<double> q2(Kl,Km);

  //--- sqrt from complex number parameter
        std::complex<double> sqrt_coef_1 = std::sqrt(costeta1+s2*sinteta1);
        std::complex<double> sqrt_coef_2 = std::sqrt(costeta1+s1*sinteta1);

        //--- expression from which the real part will be extracted
        std::complex<double> coef_x = (1.0/(s1-s2))*(s1*p2*sqrt_coef_1-s2*p1*sqrt_coef_2);
        std::complex<double> coef_y = (1.0/(s1-s2))*(s1*q2*sqrt_coef_1-s2*q1*sqrt_coef_2);

        // --- extracting real part
        double coef_x_1 = std::real(coef_x);
        double coef_y_1 = std::real(coef_y);

        // ---- computing the final displacement values
        double delx = coef * coef_x_1;
        double dely = coef * coef_y_1;
        double delz = 0.0;

        x[i][0] += delx;
        x[i][1] += dely;
        x[i][2] += delz;
      }
    }

  }
  
  // Modify ends.
  //*******************************
  
  //*******************************
  // Start of modification for gbcrack.

  if (style == GBCRACK) {
    double K1 = utils::numeric(FLERR,arg[2],false,lmp);
    double K2 = utils::numeric(FLERR,arg[3],false,lmp);
    double K3 = utils::numeric(FLERR,arg[4],false,lmp);
    double XTIP = utils::numeric(FLERR,arg[5],false,lmp);
    double YTIP = utils::numeric(FLERR,arg[6],false,lmp);
    double a = utils::numeric(FLERR,arg[7],false,lmp);
    double REG = utils::numeric(FLERR,arg[8],false,lmp);
    double P1RE = utils::numeric(FLERR,arg[9],false,lmp);
    double P1IM = utils::numeric(FLERR,arg[10],false,lmp);
    double P2RE = utils::numeric(FLERR,arg[11],false,lmp);
    double P2IM = utils::numeric(FLERR,arg[12],false,lmp);
    double P3RE = utils::numeric(FLERR,arg[13],false,lmp);
    double P3IM = utils::numeric(FLERR,arg[14],false,lmp);

    double A11RE = utils::numeric(FLERR,arg[15],false,lmp);
    double A11IM = utils::numeric(FLERR,arg[16],false,lmp);
    double A12RE = utils::numeric(FLERR,arg[17],false,lmp);
    double A12IM = utils::numeric(FLERR,arg[18],false,lmp);
    double A13RE = utils::numeric(FLERR,arg[19],false,lmp);
    double A13IM = utils::numeric(FLERR,arg[20],false,lmp);
    double A21RE = utils::numeric(FLERR,arg[21],false,lmp);
    double A21IM = utils::numeric(FLERR,arg[22],false,lmp);
    double A22RE = utils::numeric(FLERR,arg[23],false,lmp);
    double A22IM = utils::numeric(FLERR,arg[24],false,lmp);
    double A23RE = utils::numeric(FLERR,arg[25],false,lmp);
    double A23IM = utils::numeric(FLERR,arg[26],false,lmp);
    double A31RE = utils::numeric(FLERR,arg[27],false,lmp);
    double A31IM = utils::numeric(FLERR,arg[28],false,lmp);
    double A32RE = utils::numeric(FLERR,arg[29],false,lmp);
    double A32IM = utils::numeric(FLERR,arg[30],false,lmp);
    double A33RE = utils::numeric(FLERR,arg[31],false,lmp);
    double A33IM = utils::numeric(FLERR,arg[32],false,lmp);

    double B11RE = utils::numeric(FLERR,arg[33],false,lmp);
    double B11IM = utils::numeric(FLERR,arg[34],false,lmp);
    double B12RE = utils::numeric(FLERR,arg[35],false,lmp);
    double B12IM = utils::numeric(FLERR,arg[36],false,lmp);
    double B13RE = utils::numeric(FLERR,arg[37],false,lmp);
    double B13IM = utils::numeric(FLERR,arg[38],false,lmp);
    double B21RE = utils::numeric(FLERR,arg[39],false,lmp);
    double B21IM = utils::numeric(FLERR,arg[40],false,lmp);
    double B22RE = utils::numeric(FLERR,arg[41],false,lmp);
    double B22IM = utils::numeric(FLERR,arg[42],false,lmp);
    double B23RE = utils::numeric(FLERR,arg[43],false,lmp);
    double B23IM = utils::numeric(FLERR,arg[44],false,lmp);
    double B31RE = utils::numeric(FLERR,arg[45],false,lmp);
    double B31IM = utils::numeric(FLERR,arg[46],false,lmp);
    double B32RE = utils::numeric(FLERR,arg[47],false,lmp);
    double B32IM = utils::numeric(FLERR,arg[48],false,lmp);
    double B33RE = utils::numeric(FLERR,arg[49],false,lmp);
    double B33IM = utils::numeric(FLERR,arg[50],false,lmp);

    double S11 = utils::numeric(FLERR,arg[51],false,lmp);
    double S12 = utils::numeric(FLERR,arg[52],false,lmp);
    double S13 = utils::numeric(FLERR,arg[53],false,lmp);
    double S21 = utils::numeric(FLERR,arg[54],false,lmp);
    double S22 = utils::numeric(FLERR,arg[55],false,lmp);
    double S23 = utils::numeric(FLERR,arg[56],false,lmp);
    double S31 = utils::numeric(FLERR,arg[57],false,lmp);
    double S32 = utils::numeric(FLERR,arg[58],false,lmp);
    double S33 = utils::numeric(FLERR,arg[59],false,lmp);

    double PHI11 = utils::numeric(FLERR,arg[60],false,lmp);
    double PHI12 = utils::numeric(FLERR,arg[61],false,lmp);
    double PHI13 = utils::numeric(FLERR,arg[62],false,lmp);
    double PHI21 = utils::numeric(FLERR,arg[63],false,lmp);
    double PHI22 = utils::numeric(FLERR,arg[64],false,lmp);
    double PHI23 = utils::numeric(FLERR,arg[65],false,lmp);
    double PHI31 = utils::numeric(FLERR,arg[66],false,lmp);
    double PHI32 = utils::numeric(FLERR,arg[67],false,lmp);
    double PHI33 = utils::numeric(FLERR,arg[68],false,lmp);

    // some auxiliary quantities 
    const std::complex<double> imu(0.0,1.0); // imaginary unit
    const double pi = 3.1415926535897932;
    const double conv_u = 100.0; // Conversion factor from MPa*sqrt(m)*(sqrt(A)/GPa) to A

    CMAT I; 
    I.cmat[0][0] = 1.0;
    I.cmat[1][1] = 1.0;
    I.cmat[2][2] = 1.0;

    double sign;

    if (REG == 1.0) // upper grain
    {
      sign = 1.0;
    }
    else if (REG == 2.0) // lower grain
    {
      sign = -1.0;
    }
    else 
    {
      throw std::invalid_argument("Region parameter must either be 1 (upperr grain) or 2 (lower grain)!");
    }

    // assemble required quantities
    std::complex<double> p1(P1RE,P1IM);
    std::complex<double> p2(P2RE,P2IM);
    std::complex<double> p3(P3RE,P3IM);

    CMAT A1;
    A1.cmat[0][0] = {A11RE,A11IM};
    A1.cmat[0][1] = {A12RE,A12IM};
    A1.cmat[0][2] = {A13RE,A13IM};
    A1.cmat[1][0] = {A21RE,A21IM};
    A1.cmat[1][1] = {A22RE,A22IM};
    A1.cmat[1][2] = {A23RE,A23IM};
    A1.cmat[2][0] = {A31RE,A31IM};
    A1.cmat[2][1] = {A32RE,A32IM};
    A1.cmat[2][2] = {A33RE,A33IM};

    CMAT B1_inv;
    B1_inv.cmat[0][0] = {B11RE,B11IM};
    B1_inv.cmat[0][1] = {B12RE,B12IM};
    B1_inv.cmat[0][2] = {B13RE,B13IM};
    B1_inv.cmat[1][0] = {B21RE,B21IM};
    B1_inv.cmat[1][1] = {B22RE,B22IM};
    B1_inv.cmat[1][2] = {B23RE,B23IM};
    B1_inv.cmat[2][0] = {B31RE,B31IM};
    B1_inv.cmat[2][1] = {B32RE,B32IM};
    B1_inv.cmat[2][2] = {B33RE,B33IM};

    CMAT S_brev;
    S_brev.cmat[0][0] = S11;
    S_brev.cmat[0][1] = S12;
    S_brev.cmat[0][2] = S13;
    S_brev.cmat[1][0] = S21;
    S_brev.cmat[1][1] = S22;
    S_brev.cmat[1][2] = S23;
    S_brev.cmat[2][0] = S31;
    S_brev.cmat[2][1] = S32;
    S_brev.cmat[2][2] = S33;

    // input Phi_inv to spare the trouble of inverting a complex matrix
    CMAT Phi_inv;
    Phi_inv.cmat[0][0] = PHI11;
    Phi_inv.cmat[0][1] = PHI12;
    Phi_inv.cmat[0][2] = PHI13;
    Phi_inv.cmat[1][0] = PHI21;
    Phi_inv.cmat[1][1] = PHI22;
    Phi_inv.cmat[1][2] = PHI23;
    Phi_inv.cmat[2][0] = PHI31;
    Phi_inv.cmat[2][1] = PHI32;
    Phi_inv.cmat[2][2] = PHI33;

    CVEC K; // stress intensity factors
    K.cvec[0] = K2;
    K.cvec[1] = K1;
    K.cvec[2] = K3;

    CMAT S_brev2;
    S_brev2 = cMxM(S_brev, S_brev);

    CMAT S_brev3;
    S_brev3 = cMxM(S_brev2, S_brev);

    std::complex<double> tr_Sbrev2 = S_brev2.cmat[0][0] + S_brev2.cmat[1][1] + S_brev2.cmat[2][2];

    std::complex<double> beta = 1.0e-12; // to circumvent division-by-zero-infinities

    if (std::real(tr_Sbrev2) > 1.0e-12) {
      beta = sqrt(-0.5*tr_Sbrev2);
    } 

    std::complex<double> gamma = 1.0/pi*atanh(beta); 

    // Each processor loop over the nlocal atoms it owns:
      double **x = atom->x;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) { // if atoms belong to the group specified in the command

      // Coordinates
        double x1 = x[i][0] - XTIP;
        double y1 = x[i][1] - YTIP;
        double r = sqrt(x1*x1 + y1*y1);
        double teta = atan2(y1,x1);

        if (r <= 1.0e-9) // to aviod numerical issues directly at the crack tip (tolerance is arbitrary but must be low)
        {
          r = 1.0e-9;
        }

        double sinteta = sin(teta);
        double costeta = cos(teta);

        std::complex<double> rho_1 = costeta + p1*sinteta;
        std::complex<double> rho_2 = costeta + p2*sinteta;
        std::complex<double> rho_3 = costeta + p3*sinteta;

        std::complex<double> kappa_11 = 0.5*(exp(sign*gamma*pi)*std::pow((rho_1*r/(2.0*a)),(imu*gamma)) + 
                                             exp(-sign*gamma*pi)*std::pow((rho_1*r/(2.0*a)),(-imu*gamma)));

        std::complex<double> kappa_12 = 0.5*(exp(sign*gamma*pi)*std::pow((rho_2*r/(2.0*a)),(imu*gamma)) + 
                                             exp(-sign*gamma*pi)*std::pow((rho_2*r/(2.0*a)),(-imu*gamma)));

        std::complex<double> kappa_13 = 0.5*(exp(sign*gamma*pi)*std::pow((rho_3*r/(2.0*a)),(imu*gamma)) + 
                                             exp(-sign*gamma*pi)*std::pow((rho_3*r/(2.0*a)),(-imu*gamma)));

        std::complex<double> kappa_21 = 0.5*(exp(sign*gamma*pi)*std::pow((rho_1*r/(2.0*a)),(imu*gamma)) - 
                                             exp(-sign*gamma*pi)*std::pow((rho_1*r/(2.0*a)),(-imu*gamma)));

        std::complex<double> kappa_22 = 0.5*(exp(sign*gamma*pi)*std::pow((rho_2*r/(2.0*a)),(imu*gamma)) - 
                                             exp(-sign*gamma*pi)*std::pow((rho_2*r/(2.0*a)),(-imu*gamma)));

        std::complex<double> kappa_23 = 0.5*(exp(sign*gamma*pi)*std::pow((rho_3*r/(2.0*a)),(imu*gamma)) - 
                                             exp(-sign*gamma*pi)*std::pow((rho_3*r/(2.0*a)),(-imu*gamma)));

        CMAT rho_mat;
        rho_mat.cmat[0][0] = sqrt(rho_1);
        rho_mat.cmat[1][1] = sqrt(rho_2);
        rho_mat.cmat[2][2] = sqrt(rho_3);

        CMAT rho_kappa1_mat;
        rho_kappa1_mat.cmat[0][0] = kappa_11*sqrt(rho_1);
        rho_kappa1_mat.cmat[1][1] = kappa_12*sqrt(rho_2);
        rho_kappa1_mat.cmat[2][2] = kappa_13*sqrt(rho_3);

        CMAT rho_kappa2_mat;
        rho_kappa2_mat.cmat[0][0] = kappa_21*sqrt(rho_1);
        rho_kappa2_mat.cmat[1][1] = kappa_22*sqrt(rho_2);
        rho_kappa2_mat.cmat[2][2] = kappa_23*sqrt(rho_3);

        // calculating displacement vector u in a splitted manner by employing the complex linear algebra 
        // functions, defined at the beginning of the script

        CMAT u_tmp1;
        // u_tmp1 = A1 * rho_mat * B1_inv * (I + 1/beta^2 * S_brev2)
        u_tmp1 = cMxM(cMxM(cMxM(A1,rho_mat),B1_inv) , cMpM(I,(csxM(1.0/(beta*beta),S_brev2))));
        
        CMAT u_tmp2;
        // u_tmp2 = A1 * rho_kappa1_mat * B1_inv * S_brev_2
        u_tmp2 = cMxM(cMxM(cMxM(A1,rho_kappa1_mat),B1_inv),S_brev2);

        CMAT u_tmp3;
        // u_tmp3 = A1 * rho_kappa2_mat * B1_inv * S_brev_3
        u_tmp3 = cMxM(cMxM(cMxM(A1,rho_kappa2_mat),B1_inv),S_brev3);

        std::complex<double> pref = -sqrt(1.0 - beta*beta)/(beta*beta);

        CVEC u_tmp4;
        // u_tmp4 = (u_tmp1 + pref * (u_tmp2 + i/beta^2 * u_tmp3)) * Phi_inv * K
        u_tmp4 = cMxV(cMxM(cMpM(u_tmp1,csxM(pref,cMpM(u_tmp2,csxM(imu/beta,u_tmp3)))),Phi_inv),K);

        double u_x = conv_u*sqrt(2.0*r/pi)*std::real(u_tmp4.cvec[0]);
        double u_y = conv_u*sqrt(2.0*r/pi)*std::real(u_tmp4.cvec[1]);
        double u_z = conv_u*sqrt(2.0*r/pi)*std::real(u_tmp4.cvec[2]);

        // ---- computing the final displacement values
        double delx = u_x;
        double dely = u_y;
        double delz = u_z;

        x[i][0] += delx;
        x[i][1] += dely;
        x[i][2] += delz;
      }
    }
  }
  // Modify ends.
  //*******************************

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  auto irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0)
    error->warning(FLERR,"Lost atoms via displace_atoms: original {} "
                   "current {}",atom->natoms,natoms);
}

/* ----------------------------------------------------------------------
   move atoms either by specified numeric displacement or variable evaluation
------------------------------------------------------------------------- */

void DisplaceAtoms::move(int idim, char *arg, double scale)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (strstr(arg,"v_") != arg) {
    double delta = scale*utils::numeric(FLERR,arg,false,lmp);
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) x[i][idim] += delta;

  } else {
    int ivar = input->variable->find(arg+2);
    if (ivar < 0)
      error->all(FLERR,"Variable name for displace_atoms does not exist");

    if (input->variable->equalstyle(ivar)) {
      double delta = scale * input->variable->compute_equal(ivar);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) x[i][idim] += delta;
    } else if (input->variable->atomstyle(ivar)) {
      if (mvec == nullptr) memory->create(mvec,nlocal,"displace_atoms:mvec");
      input->variable->compute_atom(ivar,igroup,mvec,1,0);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) x[i][idim] += scale*mvec[i];
    } else error->all(FLERR,"Variable for displace_atoms is invalid style");
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_atoms input line
------------------------------------------------------------------------- */

void DisplaceAtoms::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal displace_atoms command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace_atoms command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace_atoms command");
  }
}
