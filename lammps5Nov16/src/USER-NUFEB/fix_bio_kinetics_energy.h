/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

