/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_ntypes.h"

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebNtypes::ComputeNufebNtypes(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ntypes command");

  vector_flag = 1;
  extvector = 0;
  int ntypes = atom->ntypes;
  size_vector = ntypes;
  memory->create(vector,ntypes,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputeNufebNtypes::~ComputeNufebNtypes()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputeNufebNtypes::compute_vector()
{
  int ntypes = atom->ntypes;
  size_vector = ntypes;
  memory->grow(vector,ntypes,"compute:vector");

  invoked_vector = update->ntimestep;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < ntypes; i++) {
    vector[i] = 0;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      vector[t-1]++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, vector, ntypes, MPI_DOUBLE, MPI_SUM, world);
}
