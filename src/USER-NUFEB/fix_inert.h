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

FixStyle(inert,FixInert)

#else

#ifndef LMP_FIX_INERT_H
#define LMP_FIX_INERT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixInert : public Fix {
 public:
	FixInert(class LAMMPS *, int, char **);
  ~FixInert();
  int setmask();
  void init();
  void post_force(int);

 private:
  char **var;
  int *ivar;
  int seed;
  class RanPark *random;
  void inert();
};

}

#endif
#endif
