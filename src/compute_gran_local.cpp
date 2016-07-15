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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_gran_local.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

enum{P1=0,P2=1,FX=2,FY=3,FZ=4,ENERGY=5,DIST=6};

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::ComputePairGranLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute gran/local command");

  int local_narg = narg - 3;

  local_flag = 1;
  nvalues = 3 + local_narg; // contact_force(0:3) + extra values

  pstyle = new int[5+local_narg];
  pindex = new int[5+local_narg];

  pstyle[0] = P1;
  pstyle[1] = P2;
  pstyle[2] = FX;
  pstyle[3] = FY;
  pstyle[4] = FZ;

  nvalues = 5; // minumum store (6): ATOM1, ATOM2, FORCE(X,Y,Z)
  for (int iarg = 0; iarg < local_narg; iarg++) {
    if (strcmp(arg[iarg+3],"dist") == 0) {
      pstyle[nvalues++] = DIST;
    } else if (strcmp(arg[iarg+3],"energy") == 0) {
      pstyle[nvalues++] = ENERGY;
    } else {
      error->all(FLERR,"PairGran/local: Invalid keyword in compute pair/local command");
    }
  }
  size_local_cols = nvalues;

  // set singleflag if need to call pair->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (pstyle[i] != DIST) singleflag = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::~ComputePairGranLocal()
{
  memory->destroy(array_local);
  memory->destroy(array);
  delete [] pstyle;
  delete [] pindex;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init()
{
  if (singleflag && force->pair == NULL) 
    error->all(FLERR,"No gran style is defined for compute gran/local");
  if (singleflag && force->pair->single_enable == 0)
    error->all(FLERR,"Gran style does not support compute gran/local");

  /*   for (int i = 0; i < nvalues; i++) */
  /*     if (pstyle[i] == PN && pindex[i] >= force->pair->single_extra) */
  /*       error->all(FLERR,"Pair style does not have extra field" */
  /* 		 " requested by compute pair/local"); */

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute pair info

  nconts = count_pairs(0);

  //nmax = nconts
  if (nconts > nmax)
    reallocate(nconts);
  size_local_rows = nconts;
  size_local_cols = nvalues;

  if(nconts)
    compute_pairs(1);
}

/* ---------------------------------------------------------------------- */

double ComputePairGranLocal::compute_scalar()
{
  invoked_local = update->ntimestep;

  // count local entries and compute pair info
  nconts = count_pairs(0);
  if (nconts > nmax) reallocate(nconts);
  size_local_rows = nconts;
  size_local_cols = nvalues;

  return 0.0;
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
   ------------------------------------------------------------------------- */

int ComputePairGranLocal::count_pairs(int flag)
{
  int i,j,m,ii,jj,inum,jnum,itype;
  double xtmp,ytmp,ztmp,delx,dely,delz,radi,radsum;
  double rsq,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  if (flag == 0) neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  m = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
     factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;
      if (newton_pair == 0 && j >= nlocal) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      if(rsq < radsum*radsum)
        m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
   ------------------------------------------------------------------------- */

int ComputePairGranLocal::compute_pairs(int flag)
{
  int i,j,m,n,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,radi,radsum;
  double rsq,energ,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  if (flag == 0) neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  Pair *pair = force->pair;

  energ = 0;
  m = 0;
  for (ii = 0; ii < inum && m < nconts; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;
      if (newton_pair == 0 && j >= nlocal) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
      if(rsq < radsum*radsum)
      {
        jtype = type[j];
        //if (singleflag)   // always true !
        energ = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
        array_local[m][0] = tag[i];
        array_local[m][1] = tag[j];
        array_local[m][2] = delx * fpair + pair->svector[0]; // fpair : normal force normalized by branch lengh
        array_local[m][3] = dely * fpair + pair->svector[1]; // svector[0:2] : tangential force
        array_local[m][4] = delz * fpair + pair->svector[2]; // del(xyz) : branch vector
        for (n = 5; n < nvalues; n++) {
          switch (pstyle[n]) {
            case DIST:
              array_local[m][n] = sqrt(rsq) - radsum;
              break;
            case ENERGY:
              array_local[m][n] = energ;
              break;
          }
        }
        pair->eng_vdwl +=energ;
        m++;
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  memory->destroy(array_local);
  memory->create(array_local,nmax,nvalues,"gran/local:array");

  size_local_rows = n;
  size_local_cols = nvalues;

}

/* ----------------------------------------------------------------------
   memory usage of local data
   ------------------------------------------------------------------------- */

double ComputePairGranLocal::memory_usage()
{
  double bytes = nmax * nvalues * sizeof(double);
  return bytes;
}
