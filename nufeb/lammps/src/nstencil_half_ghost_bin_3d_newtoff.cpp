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

#include "nstencil_half_ghost_bin_3d_newtoff.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilHalfGhostBin3dNewtoff::
NStencilHalfGhostBin3dNewtoff(LAMMPS *lmp) : NStencil(lmp)
{
  xyzflag = 1;
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilHalfGhostBin3dNewtoff::create()
{
  int i,j,k;

  nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq) {
          stencilxyz[nstencil][0] = i;
          stencilxyz[nstencil][1] = j;
          stencilxyz[nstencil][2] = k;
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
        }
}
