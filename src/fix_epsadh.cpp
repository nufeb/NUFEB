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

#include "string.h"
#include "stdlib.h"
#include "fix_epsadh.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "mpi.h"
#include "comm.h"
#include "memory.h"
#include "input.h"
#include "variable.h"


using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 10000;
/* ---------------------------------------------------------------------- */

FixEPSAdh::FixEPSAdh(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix eps adhesion command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix eps adhesion command: calling steps should be positive integer");
  flag = force->inumeric(FLERR, arg[5]);
  if (flag != 1 && flag != 2) error->all(FLERR,"Illegal fix eps adhesion command: undefined model");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var,&arg[4][2]);
}

FixEPSAdh::~FixEPSAdh()
{
  delete [] var;
}

/* ---------------------------------------------------------------------- */

int FixEPSAdh::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEPSAdh::init()
{
  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR,"Variable name for fix eps adhesion does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR,"Variable for fix eps adhesion is invalid style");

	int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
}

/* ---------------------------------------------------------------------- */

void FixEPSAdh::init_list(int id, NeighList *ptr)
{
  list = ptr;
}


void FixEPSAdh::post_force(int vflag)
{

  double ke = input->variable->compute_equal(ivar);
  int i,ii=0,j=0,jj=0,*numneigh,inum,**firstneigh;
  int jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double outerRadi,outerRadj,radsum,rsq,r, rinv;
  double ccel,ccelx ,ccely,ccelz;
  int *jlist,*ilist;
  double del;
  double **f = atom->f;
  double **x = atom->x;
  double *outerRadius = atom->outerRadius;
  double *outerMass = atom->outerMass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *mask = atom->mask;
  double epsMassi, epsMassj, massSum;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms
  
  for (ii = 0; ii < nlocal; ii++){
	  i = ilist[ii];
	  if (!(mask[i] & groupbit)) continue;
	  xtmp = x[i][0];
	  ytmp = x[i][1];
	  ztmp = x[i][2];

	  outerRadi = outerRadius[i];
	  jlist = firstneigh[i];
	  jnum = numneigh[i];
    if (type[i] == 4) {
      epsMassi = rmass[i];
    }
    else {
      epsMassi = outerMass[i];
    }

	  for (jj = 0; jj < jnum; jj++) {
		  j = jlist[jj];
		  delx = xtmp - x[j][0];
		  dely = ytmp - x[j][1];
		  delz = ztmp - x[j][2];
		  rsq = delx*delx + dely*dely + delz*delz;

		  outerRadj = outerRadius[j];
      if (type[j] == 4) {
        epsMassj = rmass[j];
      }
      else {
        epsMassj = outerMass[j];
      }
      radsum = outerRadi + outerRadj;
      massSum = epsMassi + epsMassj;

      if(flag == 1){
      	if (rsq < 4 * radsum * radsum) {
  				r = sqrt(rsq);
  				del = r - radsum;
  				rinv = 1/r;
  				ccel = -massSum*ke*del;
      	}
      }else if(flag == 2){
      	if ((rsq < 4 * radsum * ke) && (rsq > radsum)){
  				r = sqrt(rsq);
  				del = r - radsum;
  				rinv = 1/r;
  				ccel = -massSum*ke*(radsum/r)*(radsum/r);
      	}
      }else{
      	continue;
      }

  		ccelx = delx*ccel*rinv ;
  		ccely = dely*ccel*rinv ;
  		ccelz = delz*ccel*rinv ;

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
