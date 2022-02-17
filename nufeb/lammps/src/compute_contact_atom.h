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

#ifdef COMPUTE_CLASS

ComputeStyle(contact/atom,ComputeContactAtom)

#else

#ifndef LMP_COMPUTE_CONTACT_ATOM_H
#define LMP_COMPUTE_CONTACT_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeContactAtom : public Compute {
 public:
  ComputeContactAtom(class LAMMPS *, int, char **);
  ~ComputeContactAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int nmax;
  class NeighList *list;
  double *contact;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute contact/atom requires atom style sphere

Self-explanatory.

E: Compute contact/atom requires a pair style be defined

Self-explanatory.

W: More than one compute contact/atom

It is not efficient to use compute contact/atom more than once.

*/
