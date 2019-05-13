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

FixStyle(bio/nve/sphere,FixBioNVESphere)

#else

#ifndef LMP_FIX_BIO_NVE_SPHERE_H
#define LMP_FIX_BIO_NVE_SPHERE_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixBioNVESphere : public FixNVE {
 public:
  FixBioNVESphere(class LAMMPS *, int, char **);
  virtual ~FixBioNVESphere() {}
  void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  int extra;
  int dlm;

  class FixFluid *nufebFoam;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/sphere requires atom style sphere

Self-explanatory.

E: Fix nve/sphere update dipole requires atom attribute mu

An atom style with this attribute is needed.

E: Fix nve/sphere requires extended particles

This fix can only be used for particles of a finite size.
 
E: Fix nve/sphere dlm must be used with update dipole
 
The DLM algorithm can only be used in conjunction with update dipole.


*/
