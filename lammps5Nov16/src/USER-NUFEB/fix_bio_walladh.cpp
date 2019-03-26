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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "fix_bio_walladh.h"

#include <math.h>
#include <string.h>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "pointers.h"
#include "update.h"
#include "atom_vec_bio.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER};    // XYZ PLANE need to be 0,1,2
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallAhd::FixWallAhd(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 7) error->all(FLERR,"Illegal fix walladh command");

  // if (!atom->sphere_flag)
  //   error->all(FLERR,"Fix walladh requires atom style sphere");

  restart_peratom = 1;
  create_attribute = 1;

  // wall/particle coefficients

  // wallstyle args

  int n = strlen(&arg[3][2]) + 1;
  var = new char[n];
  strcpy(var,&arg[3][2]);

  int iarg = 4;
  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(FLERR,arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(FLERR,arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all(FLERR,"Illegal fix walladh command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = force->numeric(FLERR,arg[iarg+1]);
    iarg += 2;
  }
}

/* ---------------------------------------------------------------------- */

FixWallAhd::~FixWallAhd()
{
  delete [] var;
}

/* ---------------------------------------------------------------------- */

int FixWallAhd::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  // mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallAhd::init()
{

  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR,"Variable name for fix walladh does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR,"Variable for fix walladh is invalid style");


  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixWallAhd::post_force(int vflag)
{
  double kanc = input->variable->compute_equal(ivar);
  double dx,dy,dz,del1,del2,delxy,delr,rsq;

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  double wlo = lo;
  double whi = hi;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential
  // set shear if pair potential stores history

  double **x = atom->x;
  double **f = atom->f;
  double *outer_radius = avec->outer_radius;
  double *rmass = atom->rmass;
  double *outer_mass = avec->outer_mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double eps_mass;
  double delta;
  double r, rinv, ccel, ccelx, ccely, ccelz;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] == avec->eps_mask) {
      eps_mass = rmass[i];
    }
    else {
      eps_mass = outer_mass[i];
    }

    if ((mask[i] & groupbit) && eps_mass > 0) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
        del1 = x[i][0] - wlo;
        del2 = whi - x[i][0];
        if (del1 < del2) dx = del1;
        else dx = -del2;
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - wlo;
        del2 = whi - x[i][1];
        if (del1 < del2) dy = del1;
        else dy = -del2;
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - wlo;
        del2 = whi - x[i][2];
        if (del1 < del2) dz = del1;
        else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
        if (delr > outer_radius[i]) dz = cylradius;
        else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;
    //  printf("j = %i \n", i);
      //if (rsq < 4*outer_radius[i]*outer_radius[i] && rsq > 0.95*outer_radius[i]*outer_radius[i]) {
      if (rsq <= 2*outer_radius[i]*outer_radius[i]) {
        r = sqrt(rsq);
        rinv = 1.0/r;
        delta = outer_radius[i] - r;
        ccel= delta*kanc*eps_mass;
        ccelx = dx*ccel*rinv;
        ccely = dy*ccel*rinv;
        ccelz = dz*ccel*rinv;
        f[i][0] += ccelx;
        f[i][1] += ccely;
        f[i][2] += ccelz;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallAhd::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallAhd::size_restart(int nlocal)
{
  return 4;
}

/* ---------------------------------------------------------------------- */

void FixWallAhd::reset_dt()
{
  dt = update->dt;
}
