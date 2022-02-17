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

#ifdef NTOPO_CLASS

NTopoStyle(NTOPO_ANGLE_PARTIAL,NTopoAnglePartial)

#else

#ifndef LMP_TOPO_ANGLE_PARTIAL_H
#define LMP_TOPO_ANGLE_PARTIAL_H

#include "ntopo.h"

namespace LAMMPS_NS {

class NTopoAnglePartial : public NTopo {
 public:
  NTopoAnglePartial(class LAMMPS *);
  ~NTopoAnglePartial() {}
  void build();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Angle atoms %d %d %d missing on proc %d at step %ld

UNDOCUMENTED

W: Angle atoms missing at step %ld

UNDOCUMENTED

*/
