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

#ifdef ATOM_CLASS

AtomStyle(bond,AtomVecBond)

#else

#ifndef LMP_ATOM_VEC_BOND_H
#define LMP_ATOM_VEC_BOND_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecBond : public AtomVec {
 public:
  AtomVecBond(class LAMMPS *);
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  bigint memory_usage();

 private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;
  tagint *molecule;
  int **nspecial;
  tagint **special;
  int *num_bond;
  int **bond_type;
  tagint **bond_atom;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

*/
