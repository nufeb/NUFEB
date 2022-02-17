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

#include "nstencil_half_multi_2d_newton.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfMulti2dNewton::
NStencilHalfMulti2dNewton(LAMMPS *lmp) : NStencil(lmp) {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfMulti2dNewton::create()
{
  int i,j,n;
  double rsq,typesq;
  int *s;
  double *distsq;

  int ntypes = atom->ntypes;
  for (int itype = 1; itype <= ntypes; itype++) {
    typesq = cuttypesq[itype];
    s = stencil_multi[itype];
    distsq = distsq_multi[itype];
    n = 0;
    for (j = 0; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (j > 0 || (j == 0 && i > 0)) {
          rsq = bin_distance(i,j,0);
          if (rsq < typesq) {
            distsq[n] = rsq;
            s[n++] = j*mbinx + i;
          }
        }
    nstencil_multi[itype] = n;
  }
}
