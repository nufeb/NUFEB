/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_avgph.h"

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
#include "fix_bio_kinetics.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebAvgph::ComputeNufebAvgph(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute avgcon command");

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

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebAvgph::~ComputeNufebAvgph()
{
}

/* ---------------------------------------------------------------------- */

double ComputeNufebAvgph::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int global_bgrids;

  scalar = 0;

  for(int i = 0; i < kinetics->bgrids; i++){
    scalar += kinetics->sh[i];
  }

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&kinetics->bgrids, &global_bgrids, 1, MPI_INT, MPI_SUM, world);

  scalar /= (double)global_bgrids;
  scalar = -log10(scalar);
  return scalar;
}
