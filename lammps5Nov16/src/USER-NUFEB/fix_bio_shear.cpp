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

#include "fix_bio_shear.h"

#include <string.h>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixShear::FixShear(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 10) error->all(FLERR,"Illegal fix shear command");

  tmin = force->inumeric(FLERR,arg[8]);
  tmax = force->inumeric(FLERR,arg[9]);
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0 || tmin < 0 || tmax < 0) error->all(FLERR,"Illegal fix shear command: calling steps should be positive integer");
  if (tmin > tmax) error->all(FLERR, "Illegal fix shear command: maximal calling step should be greater than minimal calling step");

  if(strcmp(arg[7], "zx") == 0) dflag = 1;
  else if(strcmp(arg[7], "zy") == 0) dflag = 2;
  else error->all(FLERR,"Illegal force vector command");

  var = new char*[3];
  ivar = new int[3];

  for (int i = 0; i < 3; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }
}

FixShear::~FixShear()
{
  for (int i = 0; i < 3; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixShear::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShear::init()
{
  for (int i = 0; i < 3; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix shear does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix shear is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void FixShear::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  if (update->ntimestep < tmin) return;
  if (update->ntimestep > tmax) return;

  double **f = atom->f;
  double **x = atom->x;
  double **v = atom->v;
  double *radius = atom->radius;
  double viscosity = input->variable->compute_equal(ivar[0]);
  double shearRate = input->variable->compute_equal(ivar[1]);
  double height = input->variable->compute_equal(ivar[2]);

  //argument test:
  //printf("nevery=%i, tmin=%i, tmax=%i, viscosity=%f, rate=%f, height = %f\n", nevery, tmin, tmax, viscosity, shearRate, height);

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
  double diameter = 2 * radius[i];
  if (dflag == 1) {
    f[i][0] += MY_3PI * viscosity * diameter * (shearRate * (x[i][2] - height)-v[i][0]);
    f[i][1] += MY_3PI * viscosity * diameter * (0.0 -v[i][1]);
    f[i][2] += MY_3PI * viscosity * diameter * (0.0 -v[i][2]);

    //f[i][0] += MY_3PI * viscosity * diameter * (shearRate * (x[i][2] - height));

  } else {
    f[i][1] += MY_3PI * viscosity * diameter * (shearRate * (x[i][2] - height)-v[i][1]);
    f[i][0] += MY_3PI * viscosity * diameter * (0.0 -v[i][0]);
    f[i][2] += MY_3PI * viscosity * diameter * (0.0 -v[i][2]);


    //f[i][1] += MY_3PI * viscosity * diameter * (shearRate * (x[i][2] - height));
  }
}
}
