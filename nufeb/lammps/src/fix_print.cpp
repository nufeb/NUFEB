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

#include "fix_print.h"
#include <mpi.h>
#include <cstring>
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPrint::FixPrint(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fp(NULL), string(NULL), copy(NULL), work(NULL), var_print(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal fix print command");
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    var_print = new char[n];
    strcpy(var_print,&arg[3][2]);
    nevery = 1;
  } else {
    nevery = force->inumeric(FLERR,arg[3]);
    if (nevery <= 0) error->all(FLERR,"Illegal fix print command");
  }

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  string = new char[n];
  strcpy(string,arg[4]);

  copy = (char *) memory->smalloc(n*sizeof(char),"fix/print:copy");
  work = (char *) memory->smalloc(n*sizeof(char),"fix/print:work");
  maxcopy = maxwork = n;

  // parse optional args

  fp = NULL;
  screenflag = 1;
  char *title = NULL;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      if (me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == NULL) {
          char str[128];
          snprintf(str,128,"Cannot open fix print file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      if (strcmp(arg[iarg+1],"yes") == 0) screenflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) screenflag = 0;
      else error->all(FLERR,"Illegal fix print command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"title") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      delete [] title;
      int n = strlen(arg[iarg+1]) + 1;
      title = new char[n];
      strcpy(title,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix print command");
  }

  // print file comment line

  if (fp && me == 0) {
    if (title) fprintf(fp,"%s\n",title);
    else fprintf(fp,"# Fix print output for fix %s\n",id);
  }

  delete [] title;
}

/* ---------------------------------------------------------------------- */

FixPrint::~FixPrint()
{
  delete [] string;
  delete [] var_print;
  memory->sfree(copy);
  memory->sfree(work);

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixPrint::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrint::init()
{
  if (var_print) {
    ivar_print = input->variable->find(var_print);
    if (ivar_print < 0)
      error->all(FLERR,"Variable name for fix print timestep does not exist");
    if (!input->variable->equalstyle(ivar_print))
      error->all(FLERR,"Variable for fix print timestep is invalid style");
    next_print = static_cast<bigint>
      (input->variable->compute_equal(ivar_print));
    if (next_print <= update->ntimestep)
      error->all(FLERR,"Fix print timestep variable returned a bad timestep");
  } else {
    if (update->ntimestep % nevery)
      next_print = (update->ntimestep/nevery)*nevery + nevery;
    else
      next_print = update->ntimestep;
  }

  // add next_print to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  modify->addstep_compute_all(next_print);
}

/* ---------------------------------------------------------------------- */

void FixPrint::setup(int /* vflag */)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixPrint::end_of_step()
{
  if (update->ntimestep != next_print) return;

  // make a copy of string to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  strcpy(copy,string);
  input->substitute(copy,work,maxcopy,maxwork,0);

  if (var_print) {
    next_print = static_cast<bigint>
      (input->variable->compute_equal(ivar_print));
    if (next_print <= update->ntimestep)
      error->all(FLERR,"Fix print timestep variable returned a bad timestep");
  } else {
    next_print = (update->ntimestep/nevery)*nevery + nevery;
  }

  modify->addstep_compute(next_print);

  if (me == 0) {
    if (screenflag && screen) fprintf(screen,"%s\n",copy);
    if (screenflag && logfile) fprintf(logfile,"%s\n",copy);
    if (fp) {
      fprintf(fp,"%s\n",copy);
      fflush(fp);
    }
  }
}
