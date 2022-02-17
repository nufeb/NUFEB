/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_gas.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"
#include "force.h"

#include "fix.h"
#include "modify.h"
#include "fix_bio_kinetics.h"
#include "fix_bio_kinetics_thermo.h"
#include "bio.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebGas::ComputeNufebGas(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute gas command");

  // register fix kinetics with this class
  kinetics = NULL;
  thermo = NULL;
  nfix = modify->nfix;

  for (int i = 0; i < nfix; i++) {
    if (strcmp(modify->fix[i]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[i]);
    } else if (strcmp(modify->fix[i]->style, "kinetics/thermo") == 0) {
      thermo = static_cast<FixKineticsThermo *>(lmp->modify->fix[i]);
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"fix kinetics command is required");
  if (thermo == NULL)
    lmp->error->all(FLERR,"fix kinetics/thermo command is required");

  pressure = force->numeric(FLERR, arg[3]);

  if (pressure < 0.0)
    error->all(FLERR, "Illegal fix kinetics/thermo command: pressure");

  bio = kinetics->bio;
  vol = kinetics->stepx * kinetics->stepy * kinetics->stepz;

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebGas::~ComputeNufebGas()
{
}

/* ---------------------------------------------------------------------- */

double ComputeNufebGas::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int global_bgrids;
  double **nuR = kinetics->nur;

  scalar = 0;

  for(int i = 0; i < kinetics->bgrids; i++){
    for (int nu = 1; nu <= bio->nnu; nu++) {
      if(bio->nustate[nu] == 1 && nuR[nu][i] < 0) {
        scalar += -nuR[nu][i] * vol * 1000;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  scalar = scalar * thermo->rg * kinetics->temp / pressure;

  return scalar;
}
