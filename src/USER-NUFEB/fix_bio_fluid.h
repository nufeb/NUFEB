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

FixStyle(nufebFoam,FixFluid)

#else

#ifndef LMP_FIX_FLUID_H
#define LMP_FIX_FLUID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFluid : public Fix {
 public:
  FixFluid(class LAMMPS *, int, char **);
 ~FixFluid();
  void init();
  int setmask();
  void post_integrate();

  int dem_steps;         // # of DEM steps run in each CFD step
  int bio_steps;         // # of biological steps run in each loop
  double bio_dt;         // biological timestep
  double dem_dt;         // DEM timestep
  int nloops;            // # of loop
  int demflag;           // 0 = biological run; 1 = DEM run
  double scale_dt;       // timestep for scaling up
  int scale_nevery;      // # of steps to perform scaling up
  int scaling;           // scale flag

  double xlo, xhi, ylo, yhi, zlo, zhi, bzhi;

 private:
  char **var;
  int *ivar;

};

}

#endif
#endif


