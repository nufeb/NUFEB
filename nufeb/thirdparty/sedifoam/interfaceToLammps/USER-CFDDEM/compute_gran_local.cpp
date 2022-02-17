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

enum{DIST,ENG,FORCE,FX,FY,FZ,PN,TAG1,TAG2};

/* ---------------------------------------------------------------------- */

ComputeGranLocal::ComputeGranLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute pair/local command");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  pstyle = new int[nvalues];
  pindex = new int[nvalues];

  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"dist") == 0) pstyle[nvalues++] = DIST;
    else if (strcmp(arg[iarg],"eng") == 0) pstyle[nvalues++] = ENG;
    else if (strcmp(arg[iarg],"force") == 0) pstyle[nvalues++] = FORCE;
    else if (strcmp(arg[iarg],"fx") == 0) pstyle[nvalues++] = FX;
    else if (strcmp(arg[iarg],"fy") == 0) pstyle[nvalues++] = FY;
    else if (strcmp(arg[iarg],"fz") == 0) pstyle[nvalues++] = FZ;
    else if (strcmp(arg[iarg],"tag1") == 0) pstyle[nvalues++] = TAG1;
    else if (strcmp(arg[iarg],"tag2") == 0) pstyle[nvalues++] = TAG2;
    else if (arg[iarg][0] == 'p') {
      int n = atoi(&arg[iarg][1]);
      if (n <= 0) error->all(FLERR,
                             "Invalid keyword in compute pair/local command");
      pstyle[nvalues] = PN;
      pindex[nvalues++] = n-1;
    } else error->all(FLERR,"Invalid keyword in compute pair/local command");
  }

  // set singleflag if need to call pair->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (pstyle[i] != DIST) singleflag = 1;

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGranLocal::~ComputeGranLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
  delete [] pstyle;
  delete [] pindex;
}

/* ---------------------------------------------------------------------- */

void ComputeGranLocal::init()
{
  if (singleflag && force->pair == NULL)
    error->all(FLERR,"No pair style is defined for compute pair/local");
  if (singleflag && force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute pair/local");

  for (int i = 0; i < nvalues; i++)
    if (pstyle[i] == PN && pindex[i] >= force->pair->single_extra)
      error->all(FLERR,"Pair style does not have extra field"
                 " requested by compute pair/local");

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeGranLocal::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeGranLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute pair info

  ncount = compute_pairs(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  compute_pairs(1);
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
------------------------------------------------------------------------- */

int ComputeGranLocal::compute_pairs(int flag)
{
  int i,j,m,n,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,eng,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *ptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  if (flag == 0) neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  tagint *tag = atom->tag;

  m = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
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
      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      if (flag) {
        if (singleflag)
          eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

        if (nvalues == 1) ptr = &vector[m];
        else ptr = array[m];

        for (n = 0; n < nvalues; n++) {
          switch (pstyle[n]) {
          case DIST:
            ptr[n] = sqrt(rsq);
            break;
          case ENG:
            ptr[n] = eng;
            break;
          case FORCE:
            ptr[n] = sqrt(rsq)*fpair;
            break;
          case FX:
            ptr[n] = delx*fpair;
            break;
          case FY:
            ptr[n] = dely*fpair;
            break;
          case FZ:
            ptr[n] = delz*fpair;
            break;
          case PN:
            ptr[n] = pair->svector[pindex[n]];
            break;
          case TAG1:
            ptr[n] = tag[i];
            break;
          case TAG2:
            ptr[n] = tag[j];
            break;
          }
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeGranLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"pair/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"pair/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeGranLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
