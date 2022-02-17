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

#include "compute_displace_atom.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDisplaceAtom::ComputeDisplaceAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  displace(NULL), id_fix(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute displace/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;
  create_attribute = 1;

  // optional args

  refreshflag = 0;
  rvar = NULL;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"refresh") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute displace/atom command");
      refreshflag = 1;
      delete [] rvar;
      int n = strlen(arg[iarg+1]) + 1;
      rvar = new char[n];
      strcpy(rvar,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute displace/atom command");
  }

  // error check

  if (refreshflag) {
    ivar = input->variable->find(rvar);
    if (ivar < 0)
      error->all(FLERR,"Variable name for compute displace/atom does not exist");
    if (input->variable->atomstyle(ivar) == 0)
      error->all(FLERR,"Compute displace/atom variable "
                 "is not atom-style variable");
  }

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;

    double **x = atom->x;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
      else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  // per-atom displacement array

  nmax = nvmax = 0;
  varatom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDisplaceAtom::~ComputeDisplaceAtom()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  memory->destroy(displace);
  delete [] rvar;
  memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeDisplaceAtom::init()
{
  // set fix which stores original atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute displace/atom fix ID");
  fix = (FixStore *) modify->fix[ifix];

  if (refreshflag) {
    ivar = input->variable->find(rvar);
    if (ivar < 0)
      error->all(FLERR,"Variable name for compute displace/atom does not exist");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeDisplaceAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local displacement array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(displace);
    nmax = atom->nmax;
    memory->create(displace,nmax,4,"displace/atom:displace");
    array_atom = displace;
  }

  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int xbox,ybox,zbox;
  double dx,dy,dz;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - xoriginal[i][2];
        displace[i][0] = dx;
        displace[i][1] = dy;
        displace[i][2] = dz;
        displace[i][3] = sqrt(dx*dx + dy*dy + dz*dz);
      } else displace[i][0] = displace[i][1] =
               displace[i][2] = displace[i][3] = 0.0;

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xoriginal[i][0];
        dy = x[i][1] + h[1]*ybox + h[3]*zbox - xoriginal[i][1];
        dz = x[i][2] + h[2]*zbox - xoriginal[i][2];
        displace[i][0] = dx;
        displace[i][1] = dy;
        displace[i][2] = dz;
        displace[i][3] = sqrt(dx*dx + dy*dy + dz*dz);
      } else displace[i][0] = displace[i][1] =
               displace[i][2] = displace[i][3] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeDisplaceAtom::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

/* ----------------------------------------------------------------------
   reset per-atom storage values, based on atom-style variable evaluation
   called by dump when dump_modify refresh is set
------------------------------------------------------------------------- */

void ComputeDisplaceAtom::refresh()
{
  if (atom->nmax > nvmax) {
    nvmax = atom->nmax;
    memory->destroy(varatom);
    memory->create(varatom,nvmax,"displace/atom:varatom");
  }

  input->variable->compute_atom(ivar,igroup,varatom,1,0);

  double **xoriginal = fix->astore;
  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (varatom[i]) domain->unmap(x[i],image[i],xoriginal[i]);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDisplaceAtom::memory_usage()
{
  double bytes = nmax*4 * sizeof(double);
  bytes += nvmax * sizeof(double);
  return bytes;
}
