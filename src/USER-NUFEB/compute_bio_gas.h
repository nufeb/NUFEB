/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(gas,ComputeNufebGas)

#else

#ifndef LMP_COMPUTE_GAS_H
#define LMP_COMPUTE_GAS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeNufebGas : public Compute {
 public:
  ComputeNufebGas(class LAMMPS *, int, char **);
  ~ComputeNufebGas();
  void init() {}
  virtual double compute_scalar();

  int nfix;                  // # of Fix objects used by dump
  char **id_fix;             // their IDs
  class Fix **fix;           // list of ptrs to the Fix objects

  double vol, pressure;      // grid vol and gas pressure

  class FixKinetics *kinetics;
  class FixKineticsThermo *thermo;
  class BIO *bio;
};
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
