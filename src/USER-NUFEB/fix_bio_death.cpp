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

#include "fix_bio_death.h"

#include <string.h>

#include "atom.h"
#include "atom_vec.h"
#include "error.h"

#include "atom_vec_bio.h"
#include "force.h"
#include "input.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeath::FixDeath(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix death requires atom style bio");

  if (narg < 5)
    error->all(FLERR, "Illegal fix death command");

  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix death command");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var, &arg[4][2]);

  demflag = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (demflag != 0 && demflag != 1)
        error->all(FLERR, "Illegal fix death command: demflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix death command");
  }

  //force_reneighbor = 1;
  //next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

FixDeath::~FixDeath() {
  delete[] var;
}

/* ---------------------------------------------------------------------- */

int FixDeath::setmask() {
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeath::init() {
  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR, "Variable name for fix death does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR, "Variable for fix death is invalid style");

  dead_dia = input->variable->compute_equal(ivar);

  if (avec->type_dead == 0) {
    error->all(FLERR, "At least one initial DEAD particle is required.");
  }
}

/* ---------------------------------------------------------------------- */

void FixDeath::pre_exchange() {
  //if (next_reneighbor != update->ntimestep) return;
  if (nevery == 0)
    return;
  if (update->ntimestep % nevery)
    return;
  if (demflag)
    return;

  death();
}

/* ---------------------------------------------------------------------- */

void FixDeath::death() {
  int * const type = atom->type;
  int * const mask = atom->mask;
  double * const radius = atom->radius;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((mask[i] & groupbit) && (mask[i] != avec->mask_dead) && (mask[i] != avec->eps_mask)) {
      if (radius[i]*2 < dead_dia) {
        type[i] = avec->type_dead;
        mask[i] = avec->mask_dead;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixDeath::modify_param(int narg, char **arg) {
  if (strcmp(arg[0], "demflag") == 0) {
    if (narg != 2)
      error->all(FLERR, "Illegal fix_modify command");
    demflag = force->inumeric(FLERR, arg[1]);
    if (demflag != 1 && demflag != 0)
      error->all(FLERR, "Illegal fix_modify command: demflag");
    return 2;
  }
  return 0;
}
