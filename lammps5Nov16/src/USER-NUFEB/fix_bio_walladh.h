/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(walladh,FixWallAhd)

#else

#ifndef LMP_FIX_WALL_ADH_H
#define LMP_FIX_WALL_ADH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallAhd : public Fix {
 public:
  FixWallAhd(class LAMMPS *, int, char **);
  virtual ~FixWallAhd();
  int setmask();
  void init();
  virtual void post_force(int);

  int size_restart(int);
  int maxsize_restart();
  void reset_dt();

 protected:
  char *var;
  int ivar;
  int wallstyle,pairstyle,wiggle,wshear,axis;
  double kn,kt,gamman,gammat,xmu;
  double lo,hi,cylradius;
  double amplitude,period,omega,vshear;
  double dt;
  int nlevels_respa;
  int time_origin;

  int *touch;
  double **shear;
  int shearupdate;

  class AtomVecBio *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix wall/gran requires atom style sphere

Self-explanatory.

E: Cannot use wall in periodic dimension

Self-explanatory.

E: Cannot wiggle and shear fix wall/gran

Cannot specify both options at the same time.

E: Invalid wiggle direction for fix wall/gran

Self-explanatory.

E: Invalid shear direction for fix wall/gran

Self-explanatory.

E: Fix wall/gran is incompatible with Pair style

Must use a granular pair style to define the parameters needed for
this fix.

*/
