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

FixStyle(sedifoam,FixSedifoam)

#else

#ifndef LMP_FIX_SEDIFOAM_H
#define LMP_FIX_SEDIFOAM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSedifoam : public Fix {
 public:
  FixSedifoam(class LAMMPS *, int, char **);
 ~FixSedifoam();
  void init();
  int setmask();

  int bio_steps;         // # of biological steps run in each loop
  double bio_dt;         // biological timestep
  double dem_dt;         // DEM timestep
  int nloops;            // # of loop
  int demflag;           // 0 = biological run; 1 = DEM run


 private:
  char **var;
  int *ivar;

};

}

#endif
#endif


