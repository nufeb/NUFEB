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

#ifdef ATOM_CLASS

AtomStyle(bio,AtomVecBio)

#else

#ifndef LMP_ATOM_VEC_BIO_H
#define LMP_ATOM_VEC_BIO_H

#include "atom_vec_sphere.h"

namespace LAMMPS_NS {

class AtomVecBio : public AtomVecSphere {
 public:
  AtomVecBio(class LAMMPS *);
  ~AtomVecBio() {}
  void init();
  void grow(int);
  void create_atom(int itype, double *coord);
  void data_atom(double *, imageint, char **);
  void copy(int, int, int);
  void grow_reset();
  bigint memory_usage();

 private:
  double *sub, *o2, *nh4, *no2, *no3;
  double *virtualMass;
  double *outerRadius;
  double *outerMass;
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

E: Invalid radius in Atoms section of data file

Radius must be >= 0.0.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

*/
