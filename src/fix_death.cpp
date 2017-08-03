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

#include "fix_death.h"

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "atom_vec.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDeath::FixDeath(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix death command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix death command");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var,&arg[4][2]);

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

FixDeath::~FixDeath(){
  delete [] var;
}

/* ---------------------------------------------------------------------- */

int FixDeath::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeath::init(){
  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR,"Variable name for fix death does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR,"Variable for fix death is invalid style");
}

/* ---------------------------------------------------------------------- */

void FixDeath::pre_exchange()
{
  if (next_reneighbor != update->ntimestep) return;
  death();
}

/* ---------------------------------------------------------------------- */

void FixDeath::death()
{
  double criticalDia= input->variable->compute_equal(ivar);
  double criticalMass = 1e-18;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

//(rmass[i] < criticalMass)
  for (i = nall-1; i >= 0; i--) {
  	//delete atom
  	if((mask[i] & groupbit) && (rmass[i] < criticalMass)) {
			atom->avec->copy(nall-1,i,1);
			atom->nlocal--;
			atom->natoms--;
			continue;
  	}

    if ((mask[i] & groupbit) &&
    		(type[i] == 1 || type[i] == 2 || type[i] == 3)) {
    	if(radius[i] * 2 < criticalDia) {
    		type[i] = 6;
    		mask[i] = 65;
    	}
    }
  }

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
 // printf("nlocal=%i, all=%i \n", nlocal, atom->natoms );
  next_reneighbor += nevery;
}
