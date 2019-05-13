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
#include <stdlib.h>
#include <string.h>
#include "atom_vec_bio.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "bio.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecBio::AtomVecBio(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  radvary = 1;
  comm_x_only = 0;
  comm_f_only = 0;
  size_forward = 7;
  size_reverse = 6;
  size_border = 10;
  size_velocity = 6;
  size_data_atom = 5;
  size_data_vel = 7;
  xcol_data = 5;

  atom->sphere_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag = atom->torque_flag = 1;

  outer_mass = memory->create(outer_mass,nmax,"atom:outerMass");
  outer_radius = memory->create(outer_radius,nmax,"atom:outerRadius");;
  type_eps = 0;
  type_dead = 0;
  eps_mask = 0;
  mask_het = 0;
  mask_dead = 0;

  // instantiate BIO class
  bio = new BIO(lmp);

  for (int i = 0; i < atom->nlocal; i++) {
    outer_radius[i] = 0;
    outer_mass[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

AtomVecBio::~AtomVecBio()
{
  memory->destroy(outer_mass);
  memory->destroy(outer_radius);

  delete bio;
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::init()
{
  AtomVec::init();
  set_group_mask();

  for (int i = 0; i < atom->nlocal; i++) {
    outer_radius[i] = atom->radius[i];
    outer_mass[i] = atom->rmass[i];
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecBio::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  radius = memory->grow(atom->radius,nmax,"atom:radius");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  outer_mass = memory->grow(outer_mass,nmax,"atom:outerMass");
  outer_radius = memory->grow(outer_radius,nmax,"atom:outerRadius");
  omega = memory->grow(atom->omega,nmax,3,"atom:omega");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecBio::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  radius = atom->radius; rmass = atom->rmass;
  omega = atom->omega; torque = atom->torque;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecBio::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  radius[j] = radius[i];
  rmass[j] = rmass[i];
  outer_radius[j] = outer_radius[i];
  outer_mass[j] = outer_mass[i];
  omega[j][0] = omega[i][0];
  omega[j][1] = omega[i][1];
  omega[j][2] = omega[i][2];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_comm(int n, int *list, double *buf,
			  int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = outer_radius[j];
      buf[m++] = outer_mass[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = outer_radius[j];
      buf[m++] = outer_mass[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_comm_vel(int n, int *list, double *buf,
			      int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = outer_radius[j];
      buf[m++] = outer_mass[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
	j = list[i];
	buf[m++] = x[j][0] + dx;
	buf[m++] = x[j][1] + dy;
	buf[m++] = x[j][2] + dz;
	buf[m++] = radius[j];
	buf[m++] = rmass[j];
	buf[m++] = outer_radius[j];
	buf[m++] = outer_mass[j];
	buf[m++] = v[j][0];
	buf[m++] = v[j][1];
	buf[m++] = v[j][2];
	buf[m++] = omega[j][0];
	buf[m++] = omega[j][1];
	buf[m++] = omega[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
	j = list[i];
	buf[m++] = x[j][0] + dx;
	buf[m++] = x[j][1] + dy;
	buf[m++] = x[j][2] + dz;
	buf[m++] = radius[j];
	buf[m++] = rmass[j];
	buf[m++] = outer_radius[j];
	buf[m++] = outer_mass[j];
	if (mask[i] & deform_groupbit) {
	  buf[m++] = v[j][0] + dvx;
	  buf[m++] = v[j][1] + dvy;
	  buf[m++] = v[j][2] + dvz;
	} else {
	  buf[m++] = v[j][0];
	  buf[m++] = v[j][1];
	  buf[m++] = v[j][2];
	}
	buf[m++] = omega[j][0];
	buf[m++] = omega[j][1];
	buf[m++] = omega[j][2];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
    buf[m++] = outer_radius[j];
    buf[m++] = outer_mass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    outer_radius[i] = buf[m++];
    outer_mass[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    outer_radius[i] = buf[m++];
    outer_mass[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    outer_radius[i] = buf[m++];
    outer_mass[i] = buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_border(int n, int *list, double *buf,
			    int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = outer_radius[j];
      buf[m++] = outer_mass[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = outer_radius[j];
      buf[m++] = outer_mass[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_border_vel(int n, int *list, double *buf,
				int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = outer_radius[j];
      buf[m++] = outer_mass[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
	buf[m++] = outer_radius[j];
	buf[m++] = outer_mass[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
	buf[m++] = outer_radius[j];
	buf[m++] = outer_mass[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
    buf[m++] = outer_radius[j];
    buf[m++] = outer_mass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    outer_radius[i] = buf[m++];
    outer_mass[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}


/* ---------------------------------------------------------------------- */

void AtomVecBio::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    outer_radius[i] = buf[m++];
    outer_mass[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    outer_radius[i] = buf[m++];
    outer_mass[i] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecBio::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = outer_radius[i];
  buf[m++] = outer_mass[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBio::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  outer_radius[nlocal] = buf[m++];
  outer_mass[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecBio::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 18 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecBio::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = outer_radius[i];
  buf[m++] = outer_mass[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecBio::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  outer_radius[nlocal] = buf[m++];
  outer_mass[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecBio::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  radius[nlocal] = 0.5;
  rmass[nlocal] = 4.0*MY_PI/3.0 * radius[nlocal]*radius[nlocal]*radius[nlocal];
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  outer_radius[nlocal] = radius[nlocal];
  outer_mass[nlocal] = rmass[nlocal];

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBio::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  radius[nlocal] = 0.5 * atof(values[2]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[3]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
    rmass[nlocal] = 4.0*MY_PI/3.0 *
      radius[nlocal]*radius[nlocal]*radius[nlocal] * density;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  outer_radius[nlocal] = 0.5 * atof(values[7]);
  if (outer_radius[nlocal] < radius[nlocal]) {
    error->one(FLERR,"Outer radius must be greater than or equal to radius");
  }
  outer_mass[nlocal] = (4.0*MY_PI/3.0)*
    ((outer_radius[nlocal]*outer_radius[nlocal]*outer_radius[nlocal])
     -(radius[nlocal]*radius[nlocal]*radius[nlocal])) * 30;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecBio::data_atom_hybrid(int nlocal, char **values)
{
  radius[nlocal] = 0.5 * atof(values[0]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[1]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
    rmass[nlocal] = 4.0*MY_PI/3.0 *
      radius[nlocal]*radius[nlocal]*radius[nlocal] * density;

  return 2;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecBio::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  omega[m][0] = atof(values[3]);
  omega[m][1] = atof(values[4]);
  omega[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecBio::data_vel_hybrid(int m, char **values)
{
  omega[m][0] = atof(values[0]);
  omega[m][1] = atof(values[1]);
  omega[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBio::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = 2.0*radius[i];
    if (radius[i] == 0.0) buf[i][3] = rmass[i];
    else
      buf[i][3] = rmass[i] / (4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i]);
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecBio::pack_data_hybrid(int i, double *buf)
{
  buf[0] = 2.0*radius[i];
  if (radius[i] == 0.0) buf[1] = rmass[i];
  else buf[1] = rmass[i] / (4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i]);
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBio::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecBio::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e",buf[0],buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecBio::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = omega[i][0];
    buf[i][5] = omega[i][1];
    buf[i][6] = omega[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecBio::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = omega[i][0];
  buf[1] = omega[i][1];
  buf[2] = omega[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecBio::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecBio::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecBio::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("outerMass")) bytes += memory->usage(outer_mass,nmax);
  if (atom->memcheck("outerRadius")) bytes += memory->usage(outer_radius,nmax);
  if (atom->memcheck("omega")) bytes += memory->usage(omega,nmax,3);
  if (atom->memcheck("torque"))
    bytes += memory->usage(torque,nmax*comm->nthreads,3);

  return bytes;
}

void AtomVecBio::set_group_mask() {
  for (int i = 1; i < group->ngroup; i++) {
    if (strcmp(group->names[i],"EPS") == 0) {
      eps_mask = pow(2, i) + 1;
    } else if (strcmp(group->names[i],"HET") == 0) {
      mask_het = pow(2, i) + 1;
    } else if (strcmp(group->names[i],"DEAD") == 0) {
      mask_dead = pow(2, i) + 1;
    }
  }
}
