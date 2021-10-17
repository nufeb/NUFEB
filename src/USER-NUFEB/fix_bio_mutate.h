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

FixStyle(mutate,FixMutate)

#else

#ifndef LMP_FIX_MUTATE_H
#define LMP_FIX_MUTATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMutate : public Fix {
 public:
  FixMutate(class LAMMPS *, int, char **);
  ~FixMutate() {};
  int setmask();
  void init() {};
  void pre_exchange();
  int modify_param(int, char **);

 private:
  class AtomVecBio *avec;
  class RanPark *random;

  int demflag, imutate, tmutate, seed;
  double prob;

  void mutate();
};

}

#endif
#endif
