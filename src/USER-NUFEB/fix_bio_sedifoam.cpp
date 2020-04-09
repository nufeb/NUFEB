/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_sedifoam.h"

#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "force.h"
#include "pointers.h"
#include "update.h"
#include "input.h"
#include "variable.h"
#include "atom.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSedifoam::FixSedifoam(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix sedifoam command");

  var = new char*[3];
  ivar = new int[3];

  for (int i = 0; i < 3; i++) {
    int n = strlen(&arg[3+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3+i][2]);
  }

  bio_steps = 1;
  bio_dt = 1;
  nloops = 1;
}

FixSedifoam::~FixSedifoam()
{
  int i;
  for (i = 0; i < 3; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;
}

void FixSedifoam::init()
{

  for (int n = 0; n < 3; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix cfddem does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix cfddem is invalid style");
  }

  bio_steps = input->variable->compute_equal(ivar[0]);
  bio_dt = input->variable->compute_equal(ivar[1]);
  nloops = input->variable->compute_equal(ivar[2]);

  if (bio_steps < 0)
    error->all(FLERR, "Bio steps can not be negative");
  if (bio_dt < 0)
    error->all(FLERR, "Bio timestep can not be negative");
  if (nloops < 0)
    error->all(FLERR, "CFDDEM loop can not be negative");

  demflag = 0;
  dem_dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixSedifoam::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */


