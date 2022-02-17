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

#ifdef DIHEDRAL_CLASS

DihedralStyle(nharmonic,DihedralNHarmonic)

#else

#ifndef DIHEDRAL_NHARMONIC_H
#define DIHEDRAL_NHARMONIC_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralNHarmonic : public Dihedral {
 public:
  DihedralNHarmonic(class LAMMPS *);
  ~DihedralNHarmonic();
  void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  int *nterms;
  double **a;

  void allocate();
};

}

#endif
#endif
