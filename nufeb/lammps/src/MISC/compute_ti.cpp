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

/* ----------------------------------------------------------------------
   Contributing author: Sai Jayaraman (University of Notre Dame)
------------------------------------------------------------------------- */

#include "compute_ti.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{PAIR,TAIL,KSPACE};

/* ---------------------------------------------------------------------- */

ComputeTI::ComputeTI(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), nterms(0), which(NULL), ivar1(NULL), ivar2(NULL),
  ilo(NULL), ihi(NULL), var1(NULL), var2(NULL), pptr(NULL), pstyle(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute ti command");

  peflag = 1;
  peratom_flag = 1;
  peatomflag = 1;
  scalar_flag = 1;
  extscalar = 1;
  timeflag = 1;

  // terms come in triplets
  // changed to quadruplets to include atom type

  nterms = (narg-3) / 4;
  if (narg != 4*nterms + 3) error->all(FLERR,"Illegal compute ti command");

  which = new int[nterms];
  ivar1 = new int[nterms];
  ivar2 = new int[nterms];
  ilo = new int[nterms];
  ihi = new int[nterms];
  var1 = new char*[nterms];
  var2 = new char*[nterms];
  pptr = new Pair*[nterms];
  pstyle = new char*[nterms];

  for (int m = 0; m < nterms; m++) pstyle[m] = NULL;

  // parse keywords

  nterms = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal compute ti command");
    if (strcmp(arg[iarg],"kspace") == 0) which[nterms] = KSPACE;
    else if (strcmp(arg[iarg],"tail") == 0) which[nterms] = TAIL;
    else which[nterms] = PAIR;

    int n = strlen(arg[iarg]) + 1;
    pstyle[nterms] = new char[n];
    strcpy(pstyle[nterms],arg[iarg]);
    force->bounds(FLERR,arg[iarg+1],atom->ntypes,ilo[nterms],ihi[nterms]);
    iarg += 1;

    if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
      int n = strlen(&arg[iarg+1][2]) + 1;
      var1[nterms] = new char[n];
      strcpy(var1[nterms],&arg[iarg+1][2]);
    } else error->all(FLERR,"Illegal compute ti command");
    if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
      int n = strlen(&arg[iarg+2][2]) + 1;
      var2[nterms] = new char[n];
      strcpy(var2[nterms],&arg[iarg+2][2]);
    } else error->all(FLERR,"Illegal compute ti command");

    nterms++;
    iarg += 3;
  }
}

/* --------------------------------------------------------------------- */

ComputeTI::~ComputeTI()
{
  for (int m = 0; m < nterms; m++) {
    delete [] var1[m];
    delete [] var2[m];
    delete [] pstyle[m];
  }
  delete [] which;
  delete [] ivar1;
  delete [] ivar2;
  delete [] var1;
  delete [] var2;
  delete [] ilo;
  delete [] ihi;
  delete [] pptr;
  delete [] pstyle;
}

/* --------------------------------------------------------------------- */

void ComputeTI::init()
{
  // setup and error checks

  for (int m = 0; m < nterms; m++) {
    ivar1[m] = input->variable->find(var1[m]);
    ivar2[m] = input->variable->find(var2[m]);
    if (ivar1[m] < 0 || ivar2[m] < 0)
      error->all(FLERR,"Variable name for compute ti does not exist");
    if (!input->variable->equalstyle(ivar1[m]) ||
        !input->variable->equalstyle(ivar2[m]))
      error->all(FLERR,"Variable for compute ti is invalid style");

    if (which[m] == PAIR) {
      pptr[m] = force->pair_match(pstyle[m],1);
      if (pptr[m] == NULL)
        error->all(FLERR,"Compute ti pair style does not exist");

    } else if (which[m] == TAIL) {
      if (force->pair == NULL || force->pair->tail_flag == 0)
        error->all(FLERR,"Compute ti tail when pair style does not "
                   "compute tail corrections");

    } else if (which[m] == KSPACE) {
      if (force->kspace == NULL)
        error->all(FLERR,"Compute ti kspace style does not exist");
    }
  }
}

/* --------------------------------------------------------------------- */

double ComputeTI::compute_scalar()
{
  double eng,engall,value1,value2;

  invoked_scalar = update->ntimestep;
  if (update->eflag_global != invoked_scalar)
    error->all(FLERR,"Energy was not tallied on needed timestep");

  const int nlocal = atom->nlocal;
  const int * const type = atom->type;
  double dUdl = 0.0;

  for (int m = 0; m < nterms; m++) {
    int total_flag = 0;
    if ((ihi[m]-ilo[m])==atom->ntypes) total_flag = 1;
    eng = 0.0;
    value1 = input->variable->compute_equal(ivar1[m]);
    value2 = input->variable->compute_equal(ivar2[m]);
    if (value1 == 0.0) continue;

    if (which[m] == PAIR) {
      if (total_flag) {
        eng = pptr[m]->eng_vdwl + pptr[m]->eng_coul;
        MPI_Allreduce(&eng,&engall,1,MPI_DOUBLE,MPI_SUM,world);
      }
      else {
        int npair = nlocal;
        double *eatom = pptr[m]->eatom;

        if (force->newton_pair) npair += atom->nghost;
        for (int i = 0; i < npair; i++)
          if ((ilo[m]<=type[i])&(ihi[m]>=type[i])) eng += eatom[i];
        MPI_Allreduce(&eng,&engall,1,MPI_DOUBLE,MPI_SUM,world);
      }
      dUdl += engall/value1 * value2;

    } else if (which[m] == TAIL) {
      double vol = domain->xprd*domain->yprd*domain->zprd;
      if (total_flag)
        eng = force->pair->etail / vol;
      else {
        eng = 0;
        for (int it = 1; it <= atom->ntypes; it++) {
          int jt;
          if ((it >= ilo[m])&&(it <=ihi[m])) jt = it;
          else jt = ilo[m];
          for (; jt <=ihi[m];jt++) {
            if ((force->pair->tail_flag)&&(force->pair->setflag[it][jt])) {
              force->pair->init_one(it,jt);
              eng += force->pair->etail_ij;
            }
            if (it !=jt) eng += force->pair->etail_ij;
          }
        }
        eng /= vol;
      }
      dUdl += eng/value1 * value2;

    } else if (which[m] == KSPACE) {
      if (total_flag)
        eng = force->kspace->energy;
      else {
        double *eatom = force->kspace->eatom;
        for(int i = 0; i < nlocal; i++)
          if ((ilo[m]<=type[i])&(ihi[m]>=type[i]))
            eng += eatom[i];
        MPI_Allreduce(&eng,&engall,1,MPI_DOUBLE,MPI_SUM,world);
        eng = engall;
      }
      dUdl += eng/value1 * value2;
    }
  }

  scalar = dUdl;
  return scalar;
}
