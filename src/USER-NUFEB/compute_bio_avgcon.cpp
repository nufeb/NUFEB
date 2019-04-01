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

#include "compute_bio_avgcon.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"

#include "fix.h"
#include "modify.h"
#include "bio.h"
#include "fix_bio_kinetics.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebAvgcon::ComputeNufebAvgcon(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute avgcon command");

  vector_flag = 1;
  extvector = 0;

  // register fix kinetics with this class
  kinetics = NULL;
  nfix = modify->nfix;

  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required");

  bio = kinetics->bio;

  size_vector = bio->nnu+1;
  memory->create(vector,bio->nnu+1,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputeNufebAvgcon::~ComputeNufebAvgcon()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputeNufebAvgcon::compute_vector()
{
  invoked_vector = update->ntimestep;

  int nnus = bio->nnu;
  int global_bgrids;

  for(int i = 1; i < bio->nnu+1; i++){
    double s = 0;

    for(int j = 0; j < kinetics->bgrids; j++){
      s += kinetics->nus[i][j];
    }

    vector[i] = s;
  }

  MPI_Allreduce(MPI_IN_PLACE, vector, nnus+1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&kinetics->bgrids, &global_bgrids, 1, MPI_INT, MPI_SUM, world);

  for(int i = 1; i < bio->nnu+1; i++){
    vector[i] /= (double)global_bgrids;
  }
}
