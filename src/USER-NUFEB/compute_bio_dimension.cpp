/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_dimension.h"

#include <math.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebDimension::ComputeNufebDimension(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute fractal dimension command");

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebDimension::~ComputeNufebDimension()
{
}

/* ---------------------------------------------------------------------- */
double ComputeNufebDimension::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  scalar = 0;

  double sums[3] = {0}; // m*d, m, r

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double dia = radius[i] * 2;
      sums[0] += rmass[i] * dia * dia;
      sums[1] += rmass[i];
      sums[2] += radius[i];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_DOUBLE, MPI_SUM, world);

  double ra = pow(sums[0] / sums[1], 0.5);
  double rm = sums[2] / atom->natoms;

  scalar = log(ra / rm) / log(atom->natoms);

  return scalar;
}
