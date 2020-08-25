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
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "math_const.h"
#include "memory.h"

#include "bio.h"
#include "atom_vec_bio.h"
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

enum{HET, AOB, NOB, ANA, COM, EPS, DEAD};
/* ---------------------------------------------------------------------- */

FixKineticsMonod::FixKineticsMonod(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {

  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (narg < 3)
    error->all(FLERR, "Not enough arguments in fix kinetics/growth/monod command");

  species = NULL;
  growrate = NULL;
  
  eps_dens = 30;
  eta_het = 0;

  kinetics = NULL;

  external_gflag = 1;

  int iarg = 3;

  while (iarg < narg){
    if (strcmp(arg[iarg],"gflag") == 0) {
      external_gflag = force->inumeric(FLERR, arg[iarg+1]);
      if (external_gflag != 0 && external_gflag != 1)
        error->all(FLERR, "Illegal fix kinetics/growth/monod command: gflag");
      iarg += 2;
    } else if (strcmp(arg[iarg],"epsdens") == 0){
	eps_dens = force->numeric(FLERR, arg[iarg+1]);
	if (eps_dens <= 0)
	  error->all(FLERR, "Illegal fix kinetics/growth/monod command: eps_dens cannot be less or equal than zero");
	iarg += 2;
    } else if (strcmp(arg[iarg],"etahet") == 0){
	eta_het = force->numeric(FLERR, arg[iarg+1]);
	if (eta_het < 0)
	  error->all(FLERR, "Illegal fix kinetics/growth/monod command: eta_het cannot be less than zero");
	iarg += 2;
    } else
      error->all(FLERR, "Illegal fix kinetics/growth/monod command");
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsMonod::~FixKineticsMonod() {
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

  if (species == NULL) {
    species = memory->create(species, atom->ntypes+1, "monod:species");
    growrate = memory->create(growrate, atom->ntypes+1, 2, kinetics->ngrids, "monod:growrate");
  }

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
 initialize parameters for monod growth
 ------------------------------------------------------------------------- */
void FixKineticsMonod::init_param() {
  mu = bio->mu;
  decay = bio->decay;
  maintain = bio->maintain;
  yield = bio->yield;
  ks = bio->ks;
  ntypes = atom->ntypes;

  isub = io2 = inh4 = ino2 = ino3 = 0;

  // initialize nutrients
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
  }

  // initialize species
  for (int i = 1; i <= ntypes; i++) {
    // take the first three chars from type name;
    char *name = new char[4];
    strncpy(name, bio->tname[i], 3);
    name[3] = 0;

    if (strcmp(name, "het") == 0 || strcmp(name, "HET") == 0) {
      species[i] = HET;
      if (isub == 0) error->all(FLERR, "het growth requires nutrient 'sub' (substrate) defined in Nutrients section");
      if (io2 == 0) error->all(FLERR, "het growth requires nutrient 'o2' defined in Nutrients section");
      if (eta_het > 0) {
	if (ino2 == 0) error->all(FLERR, "het anaerobic growth requires nutrient 'no2' defined in Nutrients section");
	if (ino3 == 0) error->all(FLERR, "het anaerobic growth requires nutrient 'no3' defined in Nutrients section");
      }
    } else if (strcmp(name, "aob") == 0 || strcmp(name, "AOB") == 0) {
      species[i] = AOB;
      if (inh4 == 0) error->all(FLERR, "aob growth requires nutrient 'nh4' defined in Nutrients section");
      if (ino2 == 0) error->all(FLERR, "aob growth requires nutrient 'no2' defined in Nutrients section");
      if (io2 == 0) error->all(FLERR, "aob growth requires nutrient 'o2' defined in Nutrients section");
    } else if (strcmp(name, "nob") == 0 || strcmp(name, "NOB") == 0) {
      species[i] = NOB;
      if (ino2 == 0) error->all(FLERR, "nob growth requires nutrient 'no2' defined in Nutrients section");
      if (ino3 == 0) error->all(FLERR, "nob growth requires nutrient 'no3' defined in Nutrients section");
      if (io2 == 0) error->all(FLERR, "nob growth requires nutrient 'o2' defined in Nutrients section");
    } else if (strcmp(name, "ana") == 0 || strcmp(name, "ANA") == 0) {
      species[i] = ANA;
      if (inh4 == 0) error->all(FLERR, "anammox growth requires nutrient 'nh4' defined in Nutrients section");
      if (ino2 == 0) error->all(FLERR, "anammox growth requires nutrient 'no2' defined in Nutrients section");
      if (io2 == 0) error->all(FLERR, "anammox growth requires nutrient 'o2' defined in Nutrients section");
    } else if (strcmp(name, "com") == 0 || strcmp(name, "COM") == 0) {
      species[i] = COM;
      if (inh4 == 0) error->all(FLERR, "comammox growth requires nutrient 'nh4' defined in Nutrients section");
      if (ino3 == 0) error->all(FLERR, "comammox growth requires nutrient 'no3' defined in Nutrients section");
      if (io2 == 0) error->all(FLERR, "comammox growth requires nutrient 'o2' defined in Nutrients section");
    } else if(strcmp(name, "eps") == 0 || strcmp(name, "EPS") == 0) {
      species[i] = EPS;
      ieps = i;
    } else if(strcmp(name, "dea") == 0 || strcmp(name, "DEA") == 0) {
      species[i] = DEAD;
    } else {
      species[i] = -1;
      error->warning(FLERR, "unrecognized species found in fix_kinetics/kinetics/monod:");
      error->warning(FLERR, bio->tname[i]);
    }

    delete[] name;
  }
}

/* ---------------------------------------------------------------------- */
void FixKineticsMonod::grow_subgrid(int n) {
  growrate = memory->create(growrate, atom->ntypes+1, 2, n, "monod:growrate");
}

/* ----------------------------------------------------------------------
 metabolism and atom update
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth(double dt, int gflag) {
  radius = atom->radius;
  rmass = atom->rmass;
  outer_mass = avec->outer_mass;
  outer_radius = avec->outer_radius;

  nus = kinetics->nus;
  nur = kinetics->nur;

  xdensity = kinetics->xdensity;

  //update growth rate
  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    //grid without atom is not considered
    if(!xdensity[0][grid]) continue;

    for (int i = 1; i <= ntypes; i++) {
      growrate[i][0][grid] = 0;
      growrate[i][1][grid] = 0;
      int spec = species[i];

      if (spec == HET) {
	growth_het(i, grid);
      } else if (spec == AOB) {
	growth_aob(i, grid);
      } else if (spec == NOB) {
	growth_nob(i, grid);
      } else if (spec == ANA) {
	growth_ana(i, grid);
      } else if (spec == COM) {
	growth_com(i, grid);
      } else if (spec == EPS) {
	growth_eps(i, grid);
      } else if (spec == DEAD) {
	growth_dead(i, grid);
      }
    }
  }

  //update atom physical attributes
  if (gflag && external_gflag) update_biomass(growrate, dt);
}

/* ----------------------------------------------------------------------
 Monod growth model for heterotrophic bacteria
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_het(int i, int grid) {
  double yield_eps = 0;
  double r1, r2, r3, r4, r5, r6, r7;
  r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0;

  if (ieps != 0) yield_eps = yield[ieps];
  //het aerobic growth rate
  r1 = mu[i] * (nus[isub][grid] / (ks[i][isub] + nus[isub][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
  //het anaerobic growth rate
  if (eta_het > 0) {
    r2 = eta_het * mu[i] * (nus[isub][grid] / (ks[i][isub] + nus[isub][grid])) * (nus[ino3][grid] / (ks[i][ino3] + nus[ino3][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
    r3 = eta_het * mu[i] * (nus[isub][grid] / (ks[i][isub] + nus[isub][grid])) * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
  }
  //decay rate
  r4 = decay[i];
  //maintenance rate
  r5 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
  if (eta_het > 0) {
    r6 = (1 / 2.86) * maintain[i] * eta_het * (nus[ino3][grid] / (ks[i][ino3] + nus[ino3][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
    r7 = (1 / 1.17) * maintain[i] * eta_het * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
  }
  //nutrient utilization
  nur[isub][grid] += ((-1 / yield[i]) * ((r1 + r2 + r3) * xdensity[i][grid]));
  nur[io2][grid] += (-((1 - yield[i] - yield_eps) / yield[i]) * r1 * xdensity[i][grid]);
  nur[io2][grid] += -(r5 * xdensity[i][grid]);
  if (eta_het > 0) {
    nur[ino2][grid] += -(((1 - yield[i] - yield_eps) / (1.17 * yield[i])) * r3 * xdensity[i][grid]);
    nur[ino3][grid] += -(((1 - yield[i] - yield_eps) / (2.86 * yield[i])) * r2 * xdensity[i][grid]);
    nur[ino2][grid] += -(r7 * xdensity[i][grid]);
    nur[ino3][grid] += -(r6 * xdensity[i][grid]);
  }

  //het overall growth rate
  growrate[i][0][grid] = r1 + r2 + r3 - r4 - r5 - r6 - r7;
  //eps shell growth rate
  growrate[i][1][grid] = (yield_eps / yield[i]) * (r1 + r2 + r3);
}

/* ----------------------------------------------------------------------
 Monod growth model for ammonia-oxidizing bacteria
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_aob(int i, int grid) {
  //growth rate
  double r1 = mu[i] * (nus[inh4][grid] / (ks[i][inh4] + nus[inh4][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
  //decay rate
  double r2 = decay[i];
  //maintenance rate
  double r3 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));

  //nutrient utilization
  nur[io2][grid] += -(((3.42 - yield[i]) / yield[i]) * r1 * xdensity[i][grid]);
  // For BM3
  //nur[io2][grid] += -(((4.57 - yield[i]) / yield[i]) * R2 * xdensity[i][grid]);
  nur[inh4][grid] += -(1 / yield[i]) * r1 * xdensity[i][grid];
  nur[ino2][grid] += (1 / yield[i]) * r1 * xdensity[i][grid];
  nur[io2][grid] += -(r3 * xdensity[i][grid]);

  //aob overall growth rate
  growrate[i][0][grid] = r1 - r2 - r3;
}

/* ----------------------------------------------------------------------
 Monod growth model for nitrite-oxidizing bacteria
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_nob(int i, int grid) {
  //growth rate
  double r1 = mu[i] * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
  //decay rate
  double r2 = decay[i];
  //maintenance rate
  double r3 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));

  //nutrient utilization
  nur[io2][grid] += -(((1.15 - yield[i]) / yield[i]) * r1 * xdensity[i][grid]);
  nur[ino2][grid] += -(1 / yield[i]) * r1 * xdensity[i][grid];
  nur[ino3][grid] += (1 / yield[i]) * r1 * xdensity[i][grid];
  nur[io2][grid] += -(r3 * xdensity[i][grid]);

  //nob overall growth rate
  growrate[i][0][grid] = r1 - r2 - r3;
}

/* ----------------------------------------------------------------------
 Monod growth model for anammox bacteria
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_ana(int i, int grid) {
  //growth rate
  double r1 = mu[i] * (nus[inh4][grid] / (ks[i][inh4] + nus[inh4][grid])) * (nus[ino2][grid] / (ks[i][ino2] + nus[ino2][grid])) * (ks[i][io2] / (ks[i][io2] + nus[io2][grid]));
  //decay rate
  double r2 = decay[i];
  //maintenance rate
  double r3 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));

  //nutrient utilization
  nur[inh4][grid] +=  -(1 / yield[i]) * r1 * xdensity[i][grid];
  nur[ino2][grid] +=  -((1 / yield[i]) + (1 / 1.14)) * r1 * xdensity[i][grid];
  nur[ino3][grid] +=  (1 / 1.14) * r1 * xdensity[i][grid];

  //ana overall growth rate
  growrate[i][0][grid] = r1 - r2 -r3;
}

/* ----------------------------------------------------------------------
 Monod growth model for comammox bacteria
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_com(int i, int grid) {
  //growth rate
  double r1 = mu[i] * (nus[inh4][grid] / (ks[i][inh4] + nus[inh4][grid])) * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));
  //decay rate
  double r2 = decay[i];
  //maintenance rate
  double r3 = maintain[i] * (nus[io2][grid] / (ks[i][io2] + nus[io2][grid]));

  //nutrient utilization
  //TODO reciprocal for oxygen utilization is unknown
  nur[io2][grid] += -(((1.14 - yield[i]) / yield[i]) * r1 * xdensity[i][grid]);
  nur[inh4][grid] += -(1 / yield[i]) * r1 * xdensity[i][grid];
  nur[ino3][grid] += (1 / yield[i]) * r1 * xdensity[i][grid];
  nur[io2][grid] += -(r3 * xdensity[i][grid]);

  //ana overall growth rate
  growrate[i][0][grid] = r1 - r2 -r3;
}

/* ----------------------------------------------------------------------
 Monod growth model for extracellular polymeric substances
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_eps(int i, int grid) {
  //decay rate
  double r1 = decay[i];

  //nutrient utilization
  nur[isub][grid] += r1 * xdensity[i][grid];

  //eps overall growth rate
  growrate[i][0][grid] = -r1;
}

/* ----------------------------------------------------------------------
 Monod growth model for dead cells
 ------------------------------------------------------------------------- */
void FixKineticsMonod::growth_dead(int i, int grid) {
  //decay rate
  double r1 = decay[i];

  //nutrient utilization
  nur[isub][grid] += r1 * xdensity[i][grid];

  //eps overall growth rate
  growrate[i][0][grid] = -r1;
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixKineticsMonod::update_biomass(double ***growrate, double dt) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      int pos = kinetics->position(i);
      double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);

      rmass[i] = rmass[i] * (1 + growrate[t][0][pos] * dt);
      if (species[t] == HET && ieps != 0) {
        outer_mass[i] = four_thirds_pi * (outer_radius[i] * outer_radius[i] * outer_radius[i] - radius[i] * radius[i] * radius[i]) * eps_dens + growrate[t][1][pos] * rmass[i] * dt;
        outer_radius[i] = pow(three_quarters_pi * (rmass[i] / density + outer_mass[i] / eps_dens), third);
      } else {
        outer_mass[i] = rmass[i];
        outer_radius[i] = radius[i];
      }
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
    }
  }
}
