/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_kinetics_pressure.h"
#include "atom.h"
#include "atom_vec_bio.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "modify.h"
#include "update.h"
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixKineticsPressure::FixKineticsPressure(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (narg < 7)
    error->all(FLERR, "Not enough arguments in fix kinetics/pressure command");

  typein = force->inumeric(FLERR,arg[3]);
  groupin = group->find(arg[4]);
  if (groupin < 0)
    error->all(FLERR, "Can't find group for kinetics/pressure");
  if (strcmp(arg[5],"NULL") == 0) {
    typeout = -1;
  } else {
    typeout = force->inumeric(FLERR,arg[5]);
  }
  
  thresholdsq = force->numeric(FLERR,arg[6]);
  thresholdsq *= thresholdsq; 
  
  // create a new compute stress/atom style
  // id = fix-ID + press, compute group = all
  int n = strlen(id) + 8;
  char *id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_stress");
  char **newarg = new char*[4];
  newarg[0] = id_press;
  newarg[1] = const_cast<char *>("all");
  newarg[2] = const_cast<char *>("stress/atom");
  newarg[3] = const_cast<char *>("NULL");
  modify->add_compute(4,newarg);
  stress = modify->compute[modify->ncompute-1];
  delete [] newarg;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

FixKineticsPressure::~FixKineticsPressure() {
}

/* ---------------------------------------------------------------------- */

int FixKineticsPressure::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsPressure::init() {
  if (!atom->radius_flag)
    error->all(FLERR, "Fix requires atom attribute diameter");
}

/* ---------------------------------------------------------------------- */

void FixKineticsPressure::setup(int) {
  stress->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixKineticsPressure::initial_integrate(int) {
  stress->compute_peratom();
  double **tensor = stress->array_atom;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit || atom->mask[i] & groupin) {
      double sq = tensor[i][0] * tensor[i][0] +
	tensor[i][1] * tensor[i][1] +
	tensor[i][2] * tensor[i][2];
      double vol = 4.0 / 3.0 * MY_PI * atom->radius[i] * atom->radius[i] * atom->radius[i];
      sq /= vol * vol;
      if (sq > thresholdsq) {
	if (atom->mask[i] & groupbit) {
	  atom->type[i] = typein;
	  atom->mask[i] &= ~groupbit;
	  atom->mask[i] |= groupin;
	}
      } else {
	if (typeout > 0 && atom->type[i] == typein && atom->mask[i] & groupin) {
	  atom->type[i] = typeout;
	  atom->mask[i] &= ~groupin; 
	  atom->mask[i] |= groupbit;
	}
      }
    }
  }
  stress->addstep(update->ntimestep+1);
}
