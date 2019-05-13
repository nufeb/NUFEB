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

FixStyle(epsadh,FixEPSAdh)

#else

#ifndef LMP_FIX_EPSADH_H
#define LMP_FIX_EPSADH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEPSAdh : public Fix {

 public:
  FixEPSAdh(class LAMMPS *, int, char **);
 ~FixEPSAdh();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  virtual void post_force(int);

  double smax; //maximum seperatioin for force cutoff

 private:
  char *var;
  int ivar;
  int nmax, nvalues;
  int npairs;
  bigint laststep;
  bigint laststep_local;
  int flag;

  class NeighList *list;
  class AtomVecBio *avec;
};

}

#endif
#endif


