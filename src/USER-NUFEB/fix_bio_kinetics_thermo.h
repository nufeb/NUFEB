/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(kinetics/thermo,FixKineticsThermo)

#else

#ifndef SRC_FIX_KINETICSTHERMO_H
#define SRC_FIX_KINETICSTHERMO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsThermo : public Fix {
 public:
  FixKineticsThermo(class LAMMPS *, int, char **);
  ~FixKineticsThermo();
  void init();
  int setmask();
  void thermo(double);

  int yflag;                       // 0 = fixed yield 1 = dynamic yield
  int rflag;                       // 0 = open reactor 1 = closed reactor

  double stepx, stepy, stepz;      // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;  // simulation box size
  double vol;                      // grid volume and gas volume

  double **dgzero;
  double *khv;                     // Henry's constant
  int *liqtogas;                   // liquids convert to gas
  double gvol, rg;                 // gas volume and gas transfer constant

  class FixKinetics *kinetics;
  class BIO *bio;

  void init_dgzero();
  void init_khv();

  void dynamic_yield();
  void gas_liq_transfer(double);
  void compute_energy();
};

}

#endif
#endif

