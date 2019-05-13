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

FixStyle(kinetics/growth/energy,FixKineticsEnergy)

#else

#ifndef SRC_FIX_KINETICSENERGY_H
#define SRC_FIX_KINETICSENERGY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsEnergy : public Fix {
 public:
  FixKineticsEnergy(class LAMMPS *, int, char **);
  ~FixKineticsEnergy();
  void init();
  int setmask();
  void growth(double, int);

  double **growrate;

 private:
  char **var;
  int *ivar;

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // computational domain size
  int nx, ny, nz;                   // # of grids
  double vol;                       // grid volume

  double eps_dens;                  // EPS density
  int epsflag;                      // EPS flag

  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class BIO *bio;

 // double minimal_monod(int, int, int);
  double grid_monod(int, int);
  void update_biomass(double**, double);
};

}

#endif
#endif

