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

FixStyle(death,FixDeath)

#else

#ifndef LMP_FIX_DEATH_H
#define LMP_FIX_DEATH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeath : public Fix {
 public:
	FixDeath(class LAMMPS *, int, char **);
  ~FixDeath();
  int setmask();
  void init();
  void pre_exchange();
  int modify_param(int, char **);

 private:
  class AtomVecBio *avec;

  char *var;
  int ivar;

  int demflag;
  double dead_dia;

  void death();
};

}

#endif
#endif
