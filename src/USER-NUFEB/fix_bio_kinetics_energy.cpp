/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_kinetics_energy.h"

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
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"
#include "group.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsEnergy::FixKineticsEnergy(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (narg != 4)
    error->all(FLERR, "Not enough arguments in fix kinetics/monod command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[3 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3 + i][2]);
  }

  kinetics = NULL;
  epsflag = 0;
}

/* ---------------------------------------------------------------------- */

FixKineticsEnergy::~FixKineticsEnergy() {
  int i;
  for (i = 0; i < 1; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(growrate);
}

/* ---------------------------------------------------------------------- */

int FixKineticsEnergy::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsEnergy::init() {
  if (!atom->radius_flag)
    error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix kinetics/energy does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix kinetics/energy is invalid style");
  }

  // register fix kinetics with this class
  kinetics = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "eps_extract") == 0) {
      epsflag = 1;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for running iBM simulation");

  eps_dens = input->variable->compute_equal(ivar[0]);

  bio = kinetics->bio;

  if (bio->nnu == 0)
    error->all(FLERR, "fix_kinetics/energy requires Nutrients input");
  else if (bio->cata_coeff == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Catabolism Coeffs input");
  else if (bio->anab_coeff == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Anabolism Coeffs input");
  else if (bio->maintain == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Maintenance input");
  else if (bio->decay == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Decay input");
  else if (bio->ks == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Ks input");
  else if (bio->decay_coeff == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Decay Coeffs input");
  else if (bio->q == NULL)
    error->all(FLERR, "fix_kinetics/energy requires Consumption Rate input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  growrate = memory->create(growrate, atom->ntypes+1, kinetics->ngrids, "monod:growrate");

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
}

/* ----------------------------------------------------------------------
 metabolism and atom update
 ------------------------------------------------------------------------- */
void FixKineticsEnergy::growth(double dt, int gflag) {
  int ntypes = atom->ntypes;

  double **cata_coeff = bio->cata_coeff;
  double **anab_coeff = bio->anab_coeff;
  double *maintain = bio->maintain;
  double *decay = bio->decay;

  double **nur = kinetics->nur;
  int nnu = bio->nnu;

  double **grid_yield = kinetics->grid_yield;
  double **xdensity = kinetics->xdensity;
  int *nuconv = kinetics->nuconv;

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    //empty grid is not considered
    if(!xdensity[0][grid]) continue;

    for (int t = 1; t <= ntypes; t++) {
      double qmet, maint, inv_yield;

      qmet = bio->q[t] * grid_monod(t, grid);

      if (!kinetics->gibbs_cata[t][grid]) maint = 0;
      else maint = maintain[t] / -kinetics->gibbs_cata[t][grid];

      if (grid_yield[t][grid]) inv_yield = 1 / grid_yield[t][grid];
      else inv_yield = 0;

      for (int nu = 1; nu <= nnu; nu++) {
        if (bio->nustate[nu] != 0) continue;
        //microbe growth
        if (1.2 * maint < qmet) {
          double metCoeff = cata_coeff[t][nu] * inv_yield + anab_coeff[t][nu];
          growrate[t][grid] = grid_yield[t][grid] * (qmet - maint);
          // reaction in mol/m3
          if(!nuconv[nu]) nur[nu][grid] += growrate[t][grid] * xdensity[t][grid] * metCoeff / 24.6;
        //microbe maintenance
        } else if (qmet <= 1.2 * maint && maint <= qmet) {
          growrate[t][grid] = 0;
          if(!nuconv[nu]) nur[nu][grid] += cata_coeff[t][nu] * grid_yield[t][grid] * qmet * xdensity[t][grid] / 24.6;
        //microbe decay
        } else {
          double f;
          if (maint == 0) f = 0;
          else f = (maint - qmet) / maint;

          growrate[t][grid] = -decay[t] * f;

          if(!nuconv[nu]) nur[nu][grid] += (-growrate[t][grid] * bio->decay_coeff[t][nu] +
              cata_coeff[t][nu] * grid_yield[t][grid] * qmet) * xdensity[t][grid] / 24.6;
        }
      }
    }
  }

  if (gflag) update_biomass(growrate, dt);
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixKineticsEnergy::update_biomass(double **growrate, double dt) {
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
    int t = type[i];
    int pos = kinetics->position(i);

    double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);
    rmass[i] = rmass[i] * (1 + growrate[t][pos] * dt);


    //update mass and radius
    if (mask[i] == avec->mask_het) {
      //update HET radius
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      //update mass and radius for EPS shell if PES production is on
      //TODO
      if (epsflag == 1) {
        outer_mass[i] = four_thirds_pi * (outer_radius[i] * outer_radius[i] * outer_radius[i] - radius[i] * radius[i] * radius[i])
            * eps_dens + growrate[t][pos] * rmass[i];

        outer_radius[i] = pow(three_quarters_pi * (rmass[i] / density + outer_mass[i] / eps_dens), third);
        radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      }
    } else if (mask[i] != avec->eps_mask && mask[i] != avec->mask_dead) {
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      outer_mass[i] = rmass[i];
      outer_radius[i] = radius[i];
    }
  }
}

/* ----------------------------------------------------------------------
 get monod term w.r.t all nutrients
 ------------------------------------------------------------------------- */

double FixKineticsEnergy::grid_monod(int type, int grid) {
  double monod = 1;

  for (int i = 1; i <= bio->nnu; i++) {
    int flag = bio->ngflag[i];
    double s = kinetics->activity[i][flag][grid];
    double ks = bio->ks[type][i];

    if (ks != 0) {
      if (s <= 0) continue;
      monod *= s / (ks + s);
    }
  }

  return monod;
}

