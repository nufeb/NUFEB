/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(nugrowth,FixNuGrowth)

#else

#ifndef LMP_FIX_NUGROWTH_H
#define LMP_FIX_NUGROWTH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNuGrowth : public Fix {
 public:
  FixNuGrowth(class LAMMPS *, int, char **);
  ~FixNuGrowth();
  int setmask();
  void init();
  void pre_force(int);

 private:

  char **var;
  int *ivar;
  void change_dia();
  double compute_totalmass(int type);
  int compute_totaln(int type);
  FILE* pFile[3];
};

}

#endif
#endif
