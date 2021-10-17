/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_mutate.h"

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
#include "random_park.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMutate::FixMutate(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix mutate requires atom style bio");

  if (narg < 8)
    error->all(FLERR, "Illegal fix mutate command");

  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix mutate command");

  imutate = group->find(arg[4]);
  imutate = 1 | group->bitmask[imutate];
  tmutate = force->inumeric(FLERR, arg[5]);

  prob =  force->numeric(FLERR, arg[6]);
  seed = force->inumeric(FLERR, arg[7]);

  if (seed <= 0)
    error->all(FLERR, "Illegal fix mutate command: seed is negative");

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  demflag = 0;
  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (demflag != 0 && demflag != 1)
        error->all(FLERR, "Illegal fix mutate command: demflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix mutate command");
  }

  //force_reneighbor = 1;
  //next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

int FixMutate::setmask() {
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMutate::pre_exchange() {
  //if (next_reneighbor != update->ntimestep) return;
  if (nevery == 0)
    return;
  if (update->ntimestep % nevery)
    return;
  if (demflag)
    return;

  mutate();
}

/* ---------------------------------------------------------------------- */

void FixMutate::mutate() {
  int *mask = atom->mask;
  int *type = atom->type;
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (random->uniform() < prob) {
	     mask[i] = imutate;
	     type[i] = tmutate;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixMutate::modify_param(int narg, char **arg) {
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
