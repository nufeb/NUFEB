/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>

#include "atom.h"
#include "atom_vec_bio.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "math_const.h"
#include "memory.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "fix_bio_kinetics_monod.h"
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsMonod::FixKineticsMonod(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (narg < 5)
    error->all(FLERR, "Not enough arguments in fix kinetics/growth/monod command");

  var = new char*[2];
  ivar = new int[2];

  for (int i = 0; i < 2; i++) {
    int n = strlen(&arg[3 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3 + i][2]);
  }

  kinetics = NULL;

  external_gflag = 1;

  int iarg = 5;
  while (iarg < narg){
    if (strcmp(arg[iarg],"gflag") == 0) {
      external_gflag = force->inumeric(FLERR, arg[iarg+1]);
      if (external_gflag != 0 && external_gflag != 1)
        error->all(FLERR, "Illegal fix kinetics/growth/monod command: gflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix kinetics/growth/monod command");
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsMonod::~FixKineticsMonod() {
  int i;
  for (i = 0; i < 2; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(species);
  memory->destroy(growrate);
}

/* ---------------------------------------------------------------------- */

int FixKineticsMonod::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsMonod::init() {
  if (!atom->radius_flag)
    error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < 2; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix kinetics/monod does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix kinetics/monod is invalid style");
  }

  // register fix kinetics with this class
  kinetics = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "fix kinetics command is required for running IbM simulation");

  eps_dens = input->variable->compute_equal(ivar[0]);
  eta_het = input->variable->compute_equal(ivar[1]);

  bio = kinetics->bio;

  if (bio->nnu == 0)
    error->all(FLERR, "fix_kinetics/monod requires Nutrients input");
  else if (bio->maintain == NULL)
    error->all(FLERR, "fix_kinetics/monod requires Maintenance input");
  else if (bio->decay == NULL)
    error->all(FLERR, "fix_kinetics/monod requires Decay input");
  else if (bio->ks == NULL)
    error->all(FLERR, "fix_kinetics/monod requires Ks input");
  else if (bio->yield == NULL)
    error->all(FLERR, "fix_kinetics/monod requires Yield input");
  else if (bio->mu == NULL)
    error->all(FLERR, "fix_kinetics/monod requires Growth Rate input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species = memory->create(species, atom->ntypes+1, "monod:species");
  growrate = memory->create(growrate, atom->ntypes+1, 2, kinetics->ngrids, "monod:growrate");

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  vol = stepx * stepy * stepz;

  init_param();
}

/* ----------------------------------------------------------------------
 initialize growth parameters
 ------------------------------------------------------------------------- */
void FixKineticsMonod::init_param() {
  isub = io2 = inh4 = ino2 = ino3 = 0;
  ieps = idead = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
    if (strcmp(bio->nuname[nu], "sub") == 0)
      isub = nu;
    else if (strcmp(bio->nuname[nu], "o2") == 0)
      io2 = nu;
    else if (strcmp(bio->nuname[nu], "nh4") == 0)
      inh4 = nu;
    else if (strcmp(bio->nuname[nu], "no2") == 0)
      ino2 = nu;
    else if (strcmp(bio->nuname[nu], "no3") == 0)
      ino3 = nu;
    else
      error->all(FLERR, "unknow nutrient in fix_kinetics/kinetics/monod");
  }

  if (isub == 0)
    error->all(FLERR, "fix_kinetics/kinetics/monod requires nutrient sub (substrate)");
  if (io2 == 0)
    error->all(FLERR, "fix_kinetics/kinetics/monod requires nutrient o2");
  if (inh4 == 0)
    error->all(FLERR, "fix_kinetics/kinetics/monod requires nutrient nh4");
  if (ino2 == 0)
    error->all(FLERR, "fix_kinetics/kinetics/monod requires nutrient no2");
  if (ino3 == 0)
    error->all(FLERR, "fix_kinetics/kinetics/monod requires nutrient no3");

  // initialize type
  for (int i = 1; i <= atom->ntypes; i++) {
    if (strcmp(bio->tname[i], "eps") == 0) {
      species[i] = 4;
      ieps = i;
    } else if (strcmp(bio->tname[i], "dead") == 0) {
      species[i] = 5;
      idead = i;
    } else {
      // take first three char
      char *name = new char[4];
      strncpy(name, bio->tname[i], 3);
      name[3] = 0;

      if (strcmp(name, "het") == 0)
        species[i] = 1;
      else if (strcmp(name, "aob") == 0)
        species[i] = 2;
      else if (strcmp(name, "nob") == 0)
        species[i] = 3;
      else
        error->all(FLERR, "unknow species in fix_kinetics/kinetics/monod");

      delete[] name;
    }
  }

  if (ieps == 0)
    (error->warning(FLERR, "EPS is not defined in fix_kinetics/kinetics/monod"));
}

void FixKineticsMonod::grow_subgrid(int n) {
  growrate = memory->create(growrate, atom->ntypes + 1, 2, n, "monod:growrate");
}

/* ----------------------------------------------------------------------
 metabolism and atom update
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth(double dt, int gflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_mass = avec->outer_mass;
  double *outer_radius = avec->outer_radius;

  double *mu = bio->mu;
  double *decay = bio->decay;
  double *maintain = bio->maintain;
  double *yield = bio->yield;
  double **ks = bio->ks;

  double **nus = kinetics->nus;
  double **nur = kinetics->nur;

  double **xdensity = kinetics->xdensity;

  int *nuconv = kinetics->nuconv;
  double yield_eps = 0;

  if (ieps != 0) yield_eps = yield[ieps];

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    //empty grid is not considered
    if(!xdensity[0][grid]) continue;

    for (int i = 1; i <= ntypes; i++) {
      int spec = species[i];

      // HET monod model
      if (spec == 1) {
        double R1 = mu[i] * (nus[isub][grid] / (ks[i][isub] + nus[isub][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
        double R4 = eta_het * mu[i] * (nus[isub][grid] / (ks[i][isub] + nus[isub][grid])) * (nus[ino3][grid] / (ks[i][ino3] + nus[ino3][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
        double R5 = eta_het * mu[i] * (nus[isub][grid] / (ks[i][isub] + nus[isub][grid])) * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
        double R6 = decay[i];

        double R10 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
        double R13 = (1 / 2.86) * maintain[i] * eta_het * (nus[ino3][grid] / (ks[i][ino3] + nus[ino3][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
        double R14 = (1 / 1.17) * maintain[i] * eta_het * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));

        nur[isub][grid] += ((-1 / yield[i]) * ((R1 + R4 + R5) * xdensity[i][grid]));
      //if (xtype[i][grid] != 0) printf("nuR = %e \n", xtype[i][grid]);
        nur[io2][grid] += (-((1 - yield[i] - yield_eps) / yield[i]) * R1 * xdensity[i][grid]);
        nur[ino2][grid] += -(((1 - yield[i] - yield_eps) / (1.17 * yield[i])) * R5 * xdensity[i][grid]);
        nur[ino3][grid] += -(((1 - yield[i] - yield_eps) / (2.86 * yield[i])) * R4 * xdensity[i][grid]);
        nur[io2][grid] += -(R10 * xdensity[i][grid]);
        nur[ino2][grid] += -(R14 * xdensity[i][grid]);
        nur[ino3][grid] += -(R13 * xdensity[i][grid]);

        growrate[i][0][grid] = R1 + R4 + R5 - R6 - R10 - R13 - R14;
        growrate[i][1][grid] = (yield_eps / yield[i]) * (R1 + R4 + R5);
      } else if (spec == 2) {
        // AOB monod model
        double R2 = mu[i] * (nus[inh4][grid] / (ks[i][inh4] + nus[inh4][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
        double R7 = decay[i];
        double R11 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
        nur[io2][grid] += -(((3.42 - yield[i]) / yield[i]) * R2 * xdensity[i][grid]);
        // For BM3
        //nur[io2][grid] += -(((4.57 - yield[i]) / yield[i]) * R2 * xdensity[i][grid]);
        nur[inh4][grid] += -(1 / yield[i]) * R2 * xdensity[i][grid];
        nur[ino2][grid] += (1 / yield[i]) * R2 * xdensity[i][grid];
        nur[io2][grid] += -(R11 * xdensity[i][grid]);

        growrate[i][0][grid] = R2 - R7 - R11;
      } else if (spec == 3) {
        // NOB monod model
        double R3 = mu[i] * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
        double R8 = decay[i];
        double R12 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));

        nur[io2][grid] += -(((1.15 - yield[i]) / yield[i]) * R3 * xdensity[i][grid]);
        nur[ino2][grid] += -(1 / yield[i]) * R3 * xdensity[i][grid];
        nur[ino3][grid] += (1 / yield[i]) * R3 * xdensity[i][grid];
        nur[io2][grid] += -(R12 * xdensity[i][grid]);

        growrate[i][0][grid] = R3 - R8 - R12;
      } else if (spec == 4) {
        // EPS monod model
        double R9 = decay[i];

        nur[isub][grid] += R9 * xdensity[i][grid];
        growrate[i][0][grid] = -decay[i];
      } else if (spec == 5) {
        // DEAD monod model
        nur[isub][grid] += (decay[i] * xdensity[i][grid]);
        growrate[i][0][grid] = -decay[i];
      }
    }
  }

  if (gflag && external_gflag) update_biomass(growrate, dt);
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixKineticsMonod::update_biomass(double ***growrate, double dt) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_mass = avec->outer_mass;
  double *outer_radius = avec->outer_radius;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      int pos = kinetics->position(i);

      double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);
      rmass[i] = rmass[i] * (1 + growrate[t][0][pos] * dt);

      if (species[t] == 1) {
        outer_mass[i] = four_thirds_pi * (outer_radius[i] * outer_radius[i] * outer_radius[i] - radius[i] * radius[i] * radius[i]) * eps_dens + growrate[t][1][pos] * rmass[i] * dt;
        outer_radius[i] = pow(three_quarters_pi * (rmass[i] / density + outer_mass[i] / eps_dens), third);
        radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      } else {
        radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
        outer_mass[i] = rmass[i];
        outer_radius[i] = radius[i];
      }
    }
  }
}
