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
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include "atom_vec_electron.h"
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

static const char cite_user_eff_package[] =
  "USER-EFF package:\n\n"
  "@Article{Jaramillo-Botero11,\n"
  " author = {A. Jaramillo-Botero, J. Su, A. Qi, W. A. Goddard III},\n"
  " title = {Large-Scale, Long-Term Nonadiabatic Electron Molecular Dynamics for Describing Material Properties and Phenomena in Extreme Environments},\n"
  " journal = {J.~Comp.~Chem.},\n"
  " year =    2011,\n"
  " volume =  32,\n"
  " pages =   {497--512}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

AtomVecElectron::AtomVecElectron(LAMMPS *lmp) : AtomVec(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_eff_package);

  comm_x_only = comm_f_only = 0;

  mass_type = 1;
  molecular = 0;
  forceclearflag = 1;

  size_forward = 4;
  size_reverse = 4;
  size_border = 9;
  size_velocity = 3;
  size_data_atom = 8;
  size_data_vel = 5;
  xcol_data = 6;

  atom->ecp_flag = 0;

  atom->electron_flag = 1;
  atom->q_flag = atom->spin_flag = atom->eradius_flag =
    atom->ervel_flag = atom->erforce_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom-electron arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecElectron::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  q = memory->grow(atom->q,nmax,"atom:q");
  spin = memory->grow(atom->spin,nmax,"atom:spin");
  eradius = memory->grow(atom->eradius,nmax,"atom:eradius");
  ervel = memory->grow(atom->ervel,nmax,"atom:ervel");
  erforce = memory->grow(atom->erforce,nmax*comm->nthreads,"atom:erforce");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecElectron::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  q = atom->q;
  eradius = atom->eradius; ervel = atom->ervel; erforce = atom->erforce;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecElectron::copy(int i, int j, int delflag)
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

  q[j] = q[i];
  spin[j] = spin[i];
  eradius[j] = eradius[i];
  ervel[j] = ervel[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

void AtomVecElectron::force_clear(int n, size_t nbytes)
{
  memset(&erforce[n],0,nbytes);
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_comm(int n, int *list, double *buf,
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
      buf[m++] = eradius[j];
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
      buf[m++] = eradius[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_comm_vel(int n, int *list, double *buf,
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
      buf[m++] = eradius[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
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
        buf[m++] = eradius[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
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
        buf[m++] = eradius[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = eradius[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecElectron::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    eradius[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecElectron::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    eradius[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    eradius[i] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = erforce[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    buf[m++] = erforce[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecElectron::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    erforce[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    erforce[j] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_border(int n, int *list, double *buf,
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
      buf[m++] = q[j];
      buf[m++] = ubuf(spin[j]).d;
      buf[m++] = eradius[j];
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
      buf[m++] = q[j];
      buf[m++] = ubuf(spin[j]).d;
      buf[m++] = eradius[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_border_vel(int n, int *list, double *buf,
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
      buf[m++] = q[j];
      buf[m++] = ubuf(spin[j]).d;
      buf[m++] = eradius[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
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
    if (domain->triclinic == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = q[j];
        buf[m++] = ubuf(spin[j]).d;
        buf[m++] = eradius[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
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
        buf[m++] = q[j];
        buf[m++] = ubuf(spin[j]).d;
        buf[m++] = eradius[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = q[j];
    buf[m++] = ubuf(spin[j]).d;
    buf[m++] = eradius[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecElectron::unpack_border(int n, int first, double *buf)
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
    q[i] = buf[m++];
    spin[i] = (int) ubuf(buf[m++]).i;
    eradius[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecElectron::unpack_border_vel(int n, int first, double *buf)
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
    q[i] = buf[m++];
    spin[i] = (int) ubuf(buf[m++]).i;
    eradius[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    q[i] = buf[m++];
    spin[i] = (int) ubuf(buf[m++]).i;
    eradius[i] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecElectron::pack_exchange(int i, double *buf)
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
  buf[m++] = q[i];
  buf[m++] = ubuf(spin[i]).d;
  buf[m++] = eradius[i];
  buf[m++] = ervel[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecElectron::unpack_exchange(double *buf)
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
  q[nlocal] = buf[m++];
  spin[nlocal] = (int) ubuf(buf[m++]).i;
  eradius[nlocal] = buf[m++];
  ervel[nlocal] = buf[m++];

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

int AtomVecElectron::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 15 * nlocal;        // Associated with pack_restart

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

int AtomVecElectron::pack_restart(int i, double *buf)
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

  buf[m++] = q[i];
  buf[m++] = ubuf(spin[i]).d;
  buf[m++] = eradius[i];
  buf[m++] = ervel[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecElectron::unpack_restart(double *buf)
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

  q[nlocal] = buf[m++];
  spin[nlocal] = (int) ubuf(buf[m++]).i;
  eradius[nlocal] = buf[m++];
  ervel[nlocal] = buf[m++];

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

void AtomVecElectron::create_atom(int itype, double *coord)
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

  q[nlocal] = 0.0;
  spin[nlocal] = 1;
  eradius[nlocal] = 1.0;
  ervel[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecElectron::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;

  if (nlocal == nmax) grow(0);

  tag[nlocal] = utils::tnumeric(FLERR,values[0],true,lmp);
  type[nlocal] = utils::inumeric(FLERR,values[1],true,lmp);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  q[nlocal] = utils::numeric(FLERR,values[2],true,lmp);
  spin[nlocal] = utils::inumeric(FLERR,values[3],true,lmp);
  if (spin[nlocal] == 3) atom->ecp_flag = 1;

  eradius[nlocal] = utils::numeric(FLERR,values[4],true,lmp);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  ervel[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecElectron::data_atom_hybrid(int nlocal, char **values)
{
  q[nlocal] = utils::numeric(FLERR,values[0],true,lmp);
  spin[nlocal] = utils::inumeric(FLERR,values[1],true,lmp);
  eradius[nlocal] = utils::numeric(FLERR,values[2],true,lmp);
  if (eradius[nlocal] < 0.0)
    error->one(FLERR,"Invalid eradius in Atoms section of data file");

  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  ervel[nlocal] = 0.0;

  return 3;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecElectron::data_vel(int m, char **values)
{
  v[m][0] = utils::numeric(FLERR,values[0],true,lmp);
  v[m][1] = utils::numeric(FLERR,values[1],true,lmp);
  v[m][2] = utils::numeric(FLERR,values[2],true,lmp);
  ervel[m] = utils::numeric(FLERR,values[3],true,lmp);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecElectron::data_vel_hybrid(int m, char **values)
{
  ervel[m] = utils::numeric(FLERR,values[0],true,lmp);
  return 1;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecElectron::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = q[i];
    buf[i][3] = ubuf(spin[i]).d;
    buf[i][4] = eradius[i];
    buf[i][5] = x[i][0];
    buf[i][6] = x[i][1];
    buf[i][7] = x[i][2];
    buf[i][8] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][10] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecElectron::pack_data_hybrid(int i, double *buf)
{
  buf[0] = q[i];
  buf[1] = ubuf(spin[i]).d;
  buf[2] = eradius[i];
  return 3;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecElectron::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %-1.16e %d %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,buf[i][2],
            (int) ubuf(buf[i][3]).i,buf[i][4],buf[i][5],buf[i][6],buf[i][7],
            (int) ubuf(buf[i][8]).i,(int) ubuf(buf[i][9]).i,
            (int) ubuf(buf[i][10]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecElectron::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %d %-1.16e",buf[0],(int) ubuf(buf[1]).i,buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecElectron::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = ervel[i];
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

int AtomVecElectron::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = ervel[i];
  return 1;
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

void AtomVecElectron::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],buf[i][4]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecElectron::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e",buf[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecElectron::property_atom(char *name)
{
  if (strcmp(name,"spin") == 0) return 0;
  if (strcmp(name,"eradius") == 0) return 1;
  if (strcmp(name,"ervel") == 0) return 2;
  if (strcmp(name,"erforce") == 0) return 3;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecElectron::pack_property_atom(int index, double *buf,
                                         int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = spin[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = eradius[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = ervel[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = erforce[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecElectron::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("q")) bytes += memory->usage(q,nmax);
  if (atom->memcheck("spin")) bytes += memory->usage(spin,nmax);
  if (atom->memcheck("eradius")) bytes += memory->usage(eradius,nmax);
  if (atom->memcheck("ervel")) bytes += memory->usage(ervel,nmax);
  if (atom->memcheck("erforce"))
    bytes += memory->usage(erforce,nmax*comm->nthreads);

  return bytes;
}
