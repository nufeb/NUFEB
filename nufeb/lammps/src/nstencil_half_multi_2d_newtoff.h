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

#ifdef NSTENCIL_CLASS

NStencilStyle(half/multi/2d/newtoff,
              NStencilHalfMulti2dNewtoff,
              NS_HALF | NS_MULTI | NS_2D | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_HALF_MULTI_2D_NEWTOFF_H
#define LMP_NSTENCIL_HALF_MULTI_2D_NEWTOFF_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfMulti2dNewtoff : public NStencil {
 public:
  NStencilHalfMulti2dNewtoff(class LAMMPS *);
  ~NStencilHalfMulti2dNewtoff() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
