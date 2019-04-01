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

#include <math.h>
#include <stdlib.h> // atoi

#include "atom.h"
#include "error.h"

#include "compute_bio_rough.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"
#include "neighbor.h"
#include "pair.h"
#include "force.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebRough::ComputeNufebRough(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute average height command");

  scalar_flag = 1;
  extscalar = 0;

  nx = 0;
  ny = 0;

  // optional fields
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "nx") == 0) {
      nx = force->numeric(FLERR, arg[iarg + 1]);
      if (nx <= 0)
	error->all(FLERR,"Illegal nx value in compute roughness");
    } else if (strcmp(arg[iarg], "ny") == 0) {
      ny = force->numeric(FLERR, arg[iarg + 1]);
      if (ny <= 0)
	error->all(FLERR,"Illegal ny value in compute roughness");
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeNufebRough::init()
{
  if (nx <= 0)
    nx = static_cast<int>(domain->prd[0] / force->pair->cutforce) + 1;
  if (ny <= 0)
    ny = static_cast<int>(domain->prd[1] / force->pair->cutforce) + 1;

  stepx = domain->prd[0] / nx;
  stepy = domain->prd[1] / ny;

  // neighbor->cutneighmax is not yet initialized
  double cutneighmax = force->pair->cutforce + neighbor->skin;
  if (stepx > cutneighmax || stepy > cutneighmax)
    error->all(FLERR, "Grid step size for compute ave_height must be smaller than master list distance cutoff\n"); 
  
  grow_subgrid();
}

/* ---------------------------------------------------------------------- */

void ComputeNufebRough::grow_subgrid() {
  grid = Grid<double, 2>(Box<double, 2>(domain->boxlo, domain->boxhi), { nx, ny });
  subgrid = Subgrid<double, 2>(grid, Box<double, 2>(domain->sublo, domain->subhi));
  ReduceGrid<ComputeNufebRough>::setup();
  maxh.resize(subgrid.cell_count());
}

/* ---------------------------------------------------------------------- */

double ComputeNufebRough::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double **x = atom->x;

  std::fill(maxh.begin(), maxh.end(), 0);

  for (int i = 0; i < nlocal + nghost; i++) {
    if ((mask[i] & groupbit) && subgrid.is_inside( { x[i][0], x[i][1] }) && x[i][2] >= domain->sublo[2] && x[i][2] < domain->subhi[2]) {
      int cell = subgrid.get_index( { x[i][0], x[i][1] });
      double z = x[i][2] + atom->radius[i] - domain->boxlo[2];
      maxh[cell] = MAX(maxh[cell], z);
    }
  }

  ReduceGrid<ComputeNufebRough>::exchange();

  double ave_height = 0;
  if (domain->boxlo[2] == domain->sublo[2]) { // bottom most subdomain
    for (int i = 0; i < maxh.size(); i++) {
      ave_height += maxh[i] * stepx * stepy;
    }
  }

  ave_height = ave_height / (domain->prd[0] * domain->prd[1]);

  MPI_Allreduce(MPI_IN_PLACE, &ave_height, 1, MPI_DOUBLE, MPI_SUM, world);

  scalar = 0;
  if (domain->boxlo[2] == domain->sublo[2]) { // bottom most subdomain
    for (int i = 0; i < maxh.size(); i++) {
      scalar += ((maxh[i] - ave_height) * (maxh[i] - ave_height)) * stepx * stepy;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  scalar = scalar / (domain->prd[0] * domain->prd[1]);
  scalar = pow(scalar, 0.5);

  return scalar;
}
