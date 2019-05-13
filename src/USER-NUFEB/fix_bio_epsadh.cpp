/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_epsadh.h"

#include <math.h>
#include <string.h>

#include "atom.h"
#include "atom_vec_bio.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pointers.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 10000;
/* ---------------------------------------------------------------------- */

FixEPSAdh::FixEPSAdh(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (narg != 6)
    error->all(FLERR, "Illegal fix eps adhesion command");

  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix eps adhesion command: calling steps should be positive integer");
  flag = force->inumeric(FLERR, arg[5]);
  if (flag != 1 && flag != 2)
    error->all(FLERR, "Illegal fix eps adhesion command: undefined model");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var, &arg[4][2]);
}

FixEPSAdh::~FixEPSAdh() {
  delete[] var;
}

/* ---------------------------------------------------------------------- */

int FixEPSAdh::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEPSAdh::init() {
  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR, "Variable name for fix eps adhesion does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR, "Variable for fix eps adhesion is invalid style");

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
}

/* ---------------------------------------------------------------------- */

void FixEPSAdh::init_list(int id, NeighList *ptr) {
  list = ptr;
}

void FixEPSAdh::post_force(int vflag) {
  double ke = input->variable->compute_equal(ivar);

  double * const * const f = atom->f;
  double * const * const x = atom->x;
  double * const outerRadius = avec->outer_radius;
  double * const outerMass = avec->outer_mass;
  double * const rmass = atom->rmass;
  int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;
  int * const mask = atom->mask;
  int * const ilist = list->ilist;
  int * const numneigh = list->numneigh;
  int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (int ii = 0; ii < nlocal; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit))
      continue;
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];

    double outerRadi = outerRadius[i];
    int * const jlist = firstneigh[i];
    int jnum = numneigh[i];

    double epsMassi = 0;
    if (atom->mask[i] == avec->eps_mask)
      epsMassi = rmass[i];
    else
      epsMassi = outerMass[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;

      double outerRadj = outerRadius[j];

      double epsMassj = 0;
      if (atom->mask[j] == avec->eps_mask)
        epsMassj = rmass[j];
      else
        epsMassj = outerMass[j];

      double radsum = outerRadi + outerRadj;
      double massSum = epsMassi + epsMassj;
      double r = sqrt(rsq);
      double del = r - radsum;
      double rinv = 1 / r;

      double ccel = 0;
      if (flag == 1) {
        if (rsq < 4 * radsum * radsum && rsq > radsum * radsum) {
          ccel = -massSum * ke * del;
        }
      } else if (flag == 2) {
        if ((rsq < 4 * radsum * ke) && (rsq > radsum * radsum)) {
          ccel = -massSum * ke * (radsum / r) * (radsum / r);
        }
      }

      double ccelx = delx * ccel * rinv;
      double ccely = dely * ccel * rinv;
      double ccelz = delz * ccel * rinv;

      f[i][0] += ccelx;
      f[i][1] += ccely;
      f[i][2] += ccelz;

      if (newton_pair || j < nlocal) {
        f[j][0] -= ccelx;
        f[j][1] -= ccely;
        f[j][2] -= ccelz;
      }
    }
  }
}
