/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_diameter.h"

#include <math.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeNufebDiameter::ComputeNufebDiameter(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute floc equivalent diameter command");

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebDiameter::~ComputeNufebDiameter()
{
}

/* ---------------------------------------------------------------------- */
double ComputeNufebDiameter::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *radius = atom->radius;

  scalar = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double volume = 4 * MY_PI * radius[i] * radius[i] * radius[i] / 3;
      double x = 6 * volume / MY_PI;
      double d = pow(x, 1.0/3);
      scalar += d;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  return scalar;
}
