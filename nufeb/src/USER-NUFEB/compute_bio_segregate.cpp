/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_bio_segregate.h"

#include <cstring>

#include "atom.h"
#include "error.h"
#include "input.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebSegregate::ComputeNufebSegregate(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute segregation command");

  scalar_flag = 1;
  extscalar = 0;

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[3+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[3+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

ComputeNufebSegregate::~ComputeNufebSegregate()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

void ComputeNufebSegregate::init()
{
  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for compute segregation does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for compute segregation is invalid style");
  }

  // Initialise constants from config file
  cutoff = input->variable->compute_equal(ivar[0]);
}

/* ---------------------------------------------------------------------- */
double ComputeNufebSegregate::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  double sumPtype = 0;

  for(int i = 0; i < nlocal; i++){
    int neighbor = 0;
    int phenotype = 0;
    for(int j = 0; j < nlocal; j++){
      if(i != j){
        double xd = x[i][0] - x[j][0];
        double yd = x[i][1] - x[j][1];
        double zd = x[i][2] - x[j][2];
        double rsq = xd * xd + yd * yd + zd * zd;

        if (rsq <= cutoff * cutoff) {
          neighbor++;
          // j and i are phenotype
          if (type[i] == type[j]) phenotype++;
        }
      }
    }
    if (neighbor != 0) {
      double ratio;
      ratio = (double)phenotype / neighbor;
      sumPtype = sumPtype + ratio;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &sumPtype, 1, MPI_DOUBLE, MPI_SUM, world);

  scalar = sumPtype / atom->natoms;

  return scalar;
}
