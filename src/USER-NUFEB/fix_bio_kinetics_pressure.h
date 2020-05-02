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

FixStyle(kinetics/pressure,FixKineticsPressure)

#else

#ifndef SRC_FIX_KINETICSPRESSURE_H
#define SRC_FIX_KINETICSPRESSURE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsPressure : public Fix {
 public:
  FixKineticsPressure(class LAMMPS *, int, char **);
  ~FixKineticsPressure();
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  
 private:
  int typein, typeout;
  int groupin;
  double thresholdsq; // threshold squared

  class AtomVecBio *avec;
  class Compute *stress;
};
}

#endif
#endif

