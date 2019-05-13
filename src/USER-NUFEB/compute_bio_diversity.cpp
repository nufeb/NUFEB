/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_diversity.h"

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebDiversity::ComputeNufebDiversity(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ntypes command");

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebDiversity::~ComputeNufebDiversity()
{
}

/* ---------------------------------------------------------------------- */

double ComputeNufebDiversity::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  bigint *ntypeList = new bigint[ntypes + 1]();

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      ntypeList[t]++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, ntypeList, ntypes + 1, MPI_INT, MPI_SUM, world);

  bigint n = 0;
  for (int i = 1; i < ntypes + 1; i++) {
    n += ntypeList[i] * (ntypeList[i] - 1);
  }

  bigint m = atom->natoms * (atom->natoms - 1);
  double x = n*1.0 / m;
  scalar = 1 - x;

  delete [] ntypeList;

  return scalar;
}
