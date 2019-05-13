/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_kinetics_ph.h"

#include <math.h>
#include <string.h>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <cfloat>

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "memory.h"
#include "comm.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "modify.h"
#include "pointers.h"
#include "variable.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsPH::FixKineticsPH(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  if (narg < 4)
    error->all(FLERR, "Not enough arguments in fix kinetics/ph command");

  //set default values
  buffer_flag = 0;
  phflag = 0;
  iph = 7.0;
  phlo = 6.5;
  phhi = 9;

  if (strcmp(arg[3], "fix") == 0)
    phflag = 0;
  else if (strcmp(arg[3], "dynamic") == 0)
    phflag = 1;
  else error->all(FLERR, "Illegal ph parameter:'fix' or 'dynamic'");

  int iarg = 4;
  while (iarg < narg){
    if (strcmp(arg[iarg],"buffer") == 0) {
      buffer_flag = force->inumeric(FLERR, arg[iarg+1]);
      phlo = force->numeric(FLERR, arg[iarg + 2]);
      phhi = force->numeric(FLERR, arg[iarg + 3]);
      if (buffer_flag != 0 && buffer_flag != 1)
        error->all(FLERR, "Illegal fix kinetics/ph command: buffer_flag");
      iarg += 4;
    } else if (strcmp(arg[iarg], "ph") == 0) {
      iph = force->numeric(FLERR, arg[iarg + 1]);
      if (iph < 0.0 || iph > 14)
        error->all(FLERR, "Illegal fix kinetics/ph command: ph");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix kinetics/ph command");
  }

}

/* ---------------------------------------------------------------------- */

FixKineticsPH::~FixKineticsPH() {
  memory->destroy(keq);
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::init() {
  // register fix kinetics with this class
  kinetics = NULL;
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for kinetics/ph styles");

  bio = kinetics->bio;
  int nnus = bio->nnu;

  if (bio->nnu == 0)
    error->all(FLERR, "fix kinetics/ph requires # of Nutrients inputs");
  else if (bio->nucharge == NULL)
    error->all(FLERR, "fix_kinetics/ph requires Nutrient Charge inputs");

  keq = memory->create(keq, nnus + 1, 4, "kinetics/ph:keq");

  init_keq();
  compute_activity(0, kinetics->ngrids, iph);
}

/* ---------------------------------------------------------------------- */

int FixKineticsPH::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::solve_ph() {
  if (!phflag) compute_activity(0, kinetics->bgrids, iph);
  else dynamic_ph(0, kinetics->bgrids);
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::init_keq() {
  // water Kj/mol
  double dG0H2O = -237.18;
  int nnus = bio->nnu;
  double **nugibbs_coeff = bio->nugibbs_coeff;

  for (int i = 1; i < nnus + 1; i++) {
    for (int j = 0; j < 4; j++) {
      int k = (j == 3) ? j : j+1;

      if (nugibbs_coeff[i][k] == INF || nugibbs_coeff[i][j] == INF) {
        keq[i][j] = 0;
      } else {
        if (j > 0) keq[i][j] = exp((nugibbs_coeff[i][k] - nugibbs_coeff[i][j]) / (-kinetics->rth * kinetics->temp));
        else keq[i][0] = exp((dG0H2O + nugibbs_coeff[i][0] - nugibbs_coeff[i][1]) / (-kinetics->rth * kinetics->temp));
      }
    }
  }
}

/* ----------------------------------------------------------------------
 compute nutrient form concentration, only called when fix ph is applied
 ------------------------------------------------------------------------- */

void FixKineticsPH::compute_activity(int first, int last, double iph) {
  int nnus = bio->nnu;
  double *sh = kinetics->sh;
  double **nus = kinetics->nus;
  double ***activity = kinetics->activity;

  double *denm = memory->create(denm, nnus + 1, "kinetics:denm");
  double gSh = pow(10, -iph);
  double gSh2 = gSh * gSh;
  double gSh3 = gSh * gSh2;

  for (int k = 1; k < nnus + 1; k++) {
    denm[k] = (1 + keq[k][0]) * gSh3 + keq[k][1] * gSh2 + keq[k][2] * keq[k][3] * gSh
        + keq[k][3] * keq[k][2] * keq[k][1];
    if (denm[k] == 0) {
      lmp->error->all(FLERR, "denm returns a zero value");
    }
    double tmp[5];
    tmp[0] = keq[k][0] * gSh3 / denm[k];
    tmp[1] = gSh3 / denm[k];
    tmp[2] = gSh2 * keq[k][1] / denm[k];
    tmp[3] = gSh * keq[k][1] * keq[k][2] / denm[k];
    tmp[4] = keq[k][1] * keq[k][2] * keq[k][3] / denm[k];
    bool is_hydrogen = false;
    if (strcmp(bio->nuname[k], "h") == 0) {
      is_hydrogen = true;
    }

#pragma ivdep
#pragma vector aligned
    for (int j = first; j < last; j++) {
      sh[j] = gSh;
      // not hydrated form acitivity
      activity[k][0][j] = nus[k][j] * tmp[0];
      // fully protonated form activity
      if (is_hydrogen) {
        activity[k][1][j] = gSh;
      } else {
        activity[k][1][j] = nus[k][j] * tmp[1];
      }
      // 1st deprotonated form activity
      activity[k][2][j] = nus[k][j] * tmp[2];
      // 2nd deprotonated form activity
      activity[k][3][j] = nus[k][j] * tmp[3];
      // 3rd deprotonated form activity
      activity[k][4][j] = nus[k][j] * tmp[4];
      // if(k==1)printf("act = %e, s= %e, flag = %i \n", activity[k][1][j], nus[k][j], bio->ngflag[k]);
    }
  }
  memory->destroy(denm);
}

/* ----------------------------------------------------------------------
 buffer ph if the value is not in defined range
 ------------------------------------------------------------------------- */

void  FixKineticsPH::buffer_ph() {
  int nnus = bio->nnu;
  int ind_na, ind_cl;

  ind_na = bio->find_nuid("na");
  ind_cl = bio->find_nuid("cl");

  if (ind_na < 0 || !ind_cl < 0)
    error->all(FLERR, "buffer ph requires nutreint 'na' and 'cl'");

  int grid;
  double oldsh, evash, ph_unbuffer;
  int **nucharge = bio->nucharge;
  double ***activity = kinetics->activity;

  // always take the last grid
  grid = kinetics->ngrids - 1;
  oldsh = kinetics->sh[grid];
  // evaluate with dynamic ph
  dynamic_ph(grid, grid+1);
  evash = kinetics->sh[grid];
  ph_unbuffer = -log10(kinetics->sh[grid]);
  kinetics->sh[grid] = oldsh;

  if (ph_unbuffer < phlo || ph_unbuffer > phhi) {
    double minus = 0;
    double plus = 0;
    compute_activity(grid, grid+1, iph);

    for (int nu = 1; nu <= nnus ; nu++){
      for (int i = 0; i < 5; i++) {
        double diff, act, chr;
        act = activity[nu][i][grid];
        chr = nucharge[nu][i];

        if (chr == INF) continue;

        diff = act * chr;

        if (diff > 0) plus += diff;
        else if (diff < 0) minus -= diff;
      }
      kinetics->sh[grid] = oldsh;
    }

    kinetics->nubs[ind_na] += minus;
    kinetics->nubs[ind_cl] += plus + evash;
  }
}

/* ----------------------------------------------------------------------
 compute ph field
 ------------------------------------------------------------------------- */

inline double sum_activity(double ***activity, double **keq, double **nus, int **nucharge, double denm, double *gsh, int w, int n, int c) {
 // not hydrated form acitivity
 activity[n][0][c] = keq[n][0] / w * nus[n][c] * gsh[2] / denm;
 // fully protonated form activity
 activity[n][1][c] = nus[n][c] * gsh[2] / denm;
 // 1st deprotonated form activity
 activity[n][2][c] = nus[n][c] * gsh[1] * keq[n][1] / denm;
 // 2nd deprotonated form activity
 activity[n][3][c] = nus[n][c] * gsh[0] * keq[n][1] * keq[n][2] / denm;
 // 3rd deprotonated form activity
 activity[n][4][c] = nus[n][c] * keq[n][1] * keq[n][2] * keq[n][3] / denm;

 double tmp[5];
 tmp[0] = nucharge[n][0] * activity[n][0][c];
 tmp[1] = nucharge[n][1] * activity[n][1][c];
 tmp[2] = nucharge[n][2] * activity[n][2][c];
 tmp[3] = nucharge[n][3] * activity[n][3][c];
 tmp[4] = nucharge[n][4] * activity[n][4][c];
 return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
}

/* ---------------------------------------------------------------------- */

inline void set_gsh(double *gsh, double value) {
 gsh[0] = value;
 gsh[1] = gsh[0] * gsh[0];
 gsh[2] = gsh[1] * gsh[0];
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::dynamic_ph(int first, int last) {
  int w = 1;

  double tol = 5e-15;
  int max_iter = 100;

  int nnus = bio->nnu;

  double *fa = memory->create(fa, kinetics->ngrids, "kinetics/ph:fa");
  double *fb = memory->create(fb, kinetics->ngrids, "kinetics/ph:fb");
  double *f = memory->create(f, kinetics->ngrids, "kinetics/ph:f");
  double *df = memory->create(df, kinetics->ngrids, "kinetics/ph:df");

  double **nus = kinetics->nus;
  double temp = kinetics->temp;
  double rth = kinetics->rth;
  double ***activity = kinetics->activity;
  int **nucharge = bio->nucharge;
  double *sh = kinetics->sh;

  double a = 1e-14;
  double b = 1;

  for (int i = first; i < last; i++) {
    fa[i] = a;
    fb[i] = b;
  }

  double gsh[3];
  set_gsh(gsh, a);
  for (int k = 1; k < nnus + 1; k++) {
    double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
      + keq[k][3] * keq[k][2] * keq[k][1];
    if (denm <= 0) {
      lmp->error->all(FLERR, "denm returns a zero value");
    }
#pragma ivdep
    for (int i = first; i < last; i++) {
      fa[i] += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i);
    }
  }

  set_gsh(gsh, b);
  for (int k = 1; k < nnus + 1; k++) {
    double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
      + keq[k][3] * keq[k][2] * keq[k][1];
    if (denm <= 0) {
      lmp->error->all(FLERR, "denm returns a zero value");
    }
#pragma ivdep
    for (int i = first; i < last; i++) {
      fb[i] += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i);
    }
  }

  bool wrong = false;
  for (int i = first; i < last; i++) {
    if (fa[i] * fb[i] > 0)
      wrong = true;
  }
  if (wrong)
    lmp->error->all(FLERR, "The sum of charges returns a wrong value");

  // Newton-Raphson method
  int ipH = 1;
  while (ipH <= max_iter) {
    for (int i = first; i < last; i++) {
      f[i] = sh[i];
      df[i] = 1;
    }

    for (int k = 1; k < nnus + 1; k++) {
#pragma ivdep
      for (int i = first; i < last; i++) {
        set_gsh(gsh, sh[i]);
        double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
          + keq[k][3] * keq[k][2] * keq[k][1];
        f[i] += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i);

        double ddenm = denm * denm;
        double aux = 3 * gsh[1] * (keq[k][0] / w + 1) + 2 * gsh[0] * keq[k][1] + keq[k][1] * keq[k][2];
        double tmp[5];
        tmp[0] = nucharge[k][0] * ((3 * gsh[1] * keq[k][0] * nus[k][i]) / (w * denm) - (keq[k][0] * nus[k][i] * gsh[2] * aux) / (w * ddenm));
        tmp[1] = nucharge[k][1] * ((3 * gsh[1] * nus[k][i]) / denm - (nus[k][i] * gsh[2] * aux) / ddenm);
        tmp[2] = nucharge[k][2] * ((2 * gsh[0] * keq[k][1] * nus[k][i]) / denm - (keq[k][1] * nus[k][i] * gsh[1] * aux) / ddenm);
        tmp[3] = nucharge[k][3] * ((keq[k][1] * keq[k][2] * nus[k][i]) / denm - (keq[k][1] * keq[k][2] * nus[k][i] * gsh[0] * aux) / ddenm);
        tmp[4] = nucharge[k][4] * (-(keq[k][1] * keq[k][2] * keq[k][3] * nus[k][i] * aux) / ddenm);
        df[i] += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
      }
    }

    // Check for convergence
    bool flag = true;
    for (int i = first; i < last; i++) {
      if (fabs(f[i]) >= tol)
      {
        flag = false;
      }
    }
    if (flag) break;

    // Compute next value
#pragma ivdep
    for (int i = first; i < last; i++) {
      if (fabs(f[i]) >= tol) {
        double d = f[i] / df[i];
        // Prevent sh below 1e-14. That can happen because sometimes the Newton
        // method overshoots to a negative sh value, due to a small derivative
        // value.
        if (d >= sh[i] - 1e-14)
          d = sh[i] / 2;
        sh[i] -= d;
      }
    }

    ipH++;
  }

  int id = bio->find_nuid("h");

  if (id > 0) {
    for (int i = first; i < last; i++) {
      activity[id][1][i] = sh[i];
    }
  }

  memory->destroy(fa);
  memory->destroy(fb);
  memory->destroy(f);
  memory->destroy(df);
}

