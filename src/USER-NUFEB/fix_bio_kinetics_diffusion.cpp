/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_kinetics_diffusion.h"

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "atom.h"
#include "domain.h"
#include "error.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "atom_vec_bio.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "memory.h"
#include "update.h"
#include "variable.h"
#include "group.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

enum{MOL, KG};
enum{PP, DD, ND, NN, DN};
enum{REGULAR, BOUNDARY, GHOST};

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion::FixKineticsDiffusion(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {

  if (narg < 8)
    error->all(FLERR, "Not enough arguments in fix diffusion command");

  shearflag = dragflag = 0;
  bulkflag = 0;
  srate = 0;
  dcflag = 0;

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[3 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3 + i][2]);
  }

  //set boundary condition flag:
  //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  if (strcmp(arg[4], "pp") == 0)
    xbcflag = PP;
  else if (strcmp(arg[4], "dd") == 0)
    xbcflag = DD;
  else if (strcmp(arg[4], "nd") == 0)
    xbcflag = ND;
  else if (strcmp(arg[4], "nn") == 0)
    xbcflag = NN;
  else if (strcmp(arg[4], "dn") == 0)
    xbcflag = DN;
  else
    error->all(FLERR, "Illegal x-axis boundary condition command");

  if (strcmp(arg[5], "pp") == 0)
    ybcflag = PP;
  else if (strcmp(arg[5], "dd") == 0)
    ybcflag = DD;
  else if (strcmp(arg[5], "nd") == 0)
    ybcflag = ND;
  else if (strcmp(arg[5], "nn") == 0)
    ybcflag = NN;
  else if (strcmp(arg[5], "dn") == 0)
    ybcflag = DN;
  else
    error->all(FLERR, "Illegal y-axis boundary condition command");

  if (strcmp(arg[6], "pp") == 0)
    zbcflag = PP;
  else if (strcmp(arg[6], "dd") == 0)
    zbcflag = DD;
  else if (strcmp(arg[6], "nd") == 0)
    zbcflag = ND;
  else if (strcmp(arg[6], "nn") == 0)
    zbcflag = NN;
  else if (strcmp(arg[6], "dn") == 0)
    zbcflag = DN;
  else
    error->all(FLERR, "Illegal z-axis boundary condition command");

  if (strcmp(arg[7], "kg") == 0)
    unit = KG;
  else if (strcmp(arg[7], "mol") == 0)
    unit = MOL;
  else
    error->all(FLERR, "Illegal unit in fix kinetics/diffusion command: specify 'kg' or 'mol'");

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "srate") == 0) {
      srate = force->numeric(FLERR, arg[iarg + 1]);
      if (srate < 0)
        error->all(FLERR, "Illegal fix kinetics/diffusion command: srate");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dcflag") == 0) {
      dcflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (dcflag < 0 || dcflag > 2)
        error->all(FLERR, "Illegal fix kinetics/diffusion command: dcflag");
      iarg += 2;
    } else if (strcmp(arg[iarg], "bulk") == 0) {
      bulkflag = 1;
      q = force->numeric(FLERR, arg[iarg + 1]);
      if (q < 0)
        error->all(FLERR, "Flow rate (frate) cannot be negative");
      rvol = force->numeric(FLERR, arg[iarg + 2]);
      if (rvol <= 0)
        error->all(FLERR, "Reactor volume (rvol) cannot be equal or less than 0");
      af = force->numeric(FLERR, arg[iarg + 3]);
      if (af < 0)
        lmp->error->all(FLERR, "Biofilm surface area (Af) cannot be negative");
      iarg += 4;
    } else
      error->all(FLERR, "Illegal fix kinetics/diffusion command");
  }
  
  setup_exchange_flag = false;
}

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion::~FixKineticsDiffusion() {
  for (int i = 0; i < 1; i++) {
    delete[] var[i];
  }

  delete[] var;
  delete[] ivar;

  memory->destroy(xgrid);
  memory->destroy(nugrid);
  memory->destroy(nuprev);
  memory->destroy(grid_diff_coeff);
  memory->destroy(ghost);

  delete[] requests;
}

/* ---------------------------------------------------------------------- */

int FixKineticsDiffusion::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsDiffusion::init() {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (!atom->radius_flag)
    error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix diffusion does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix diffusion is invalid style");
  }

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "shear") == 0) {
      shearflag = 1;
    } else if (strcmp(modify->fix[j]->style, "fdrag") == 0) {
      dragflag = 1;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for running iBM simulation");

  bio = kinetics->bio;
  tol = input->variable->compute_equal(ivar[0]);

  //set diffusion grid size
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  vol = stepx * stepy * stepz;
  bzhi = kinetics->bnz * stepz;

  if (!is_equal(stepx, stepy, stepz))
    error->all(FLERR, "Grid is not cubic");

  snxyz = kinetics->subn[0] * kinetics->subn[1] * kinetics->subn[2];

  snxx = kinetics->subn[0] + 2;
  snyy = kinetics->subn[1] + 2;
  snzz = kinetics->subn[2] + 2;

  snxx_yy = snxx * snyy;
  snxx_yy_zz = snxx_yy * snzz;
  
  int nnus = bio->nnu;

  //inlet concentration and maximum boundary condition conc value

  xgrid = memory->create(xgrid, snxx_yy_zz, 3, "diffusion:xgrid");
  nugrid = memory->create(nugrid, nnus + 1, snxx_yy_zz, "diffusion:nugrid");
  nuprev = memory->create(nuprev, nnus + 1, snxx_yy_zz, "diffusion:nuprev");
  ghost = memory->create(ghost, snxx_yy_zz, "diffusion:ghost");
  grid_diff_coeff = memory->create(grid_diff_coeff, nnus + 1, snxx_yy_zz, "diffusion:grid_diff_coeff");

  init_grid();

  // create request vector
  requests = new MPI_Request[MAX(2 * comm->nprocs, nnus + 1)];

  setup_exchange(kinetics->grid, kinetics->subgrid.get_box(), { xbcflag == 0, ybcflag == 0, zbcflag == 0 });
}

/* ----------------------------------------------------------------------
 solve diffusion and reaction
 ------------------------------------------------------------------------- */

int *FixKineticsDiffusion::diffusion(int *nuConv, int iter, double diff_dt) {
  int nnus = bio->nnu;
  this->diff_dt = diff_dt;
  double **nur = kinetics->nur;
  double **nus = kinetics->nus;
  double *nubs = kinetics->nubs;
  double **ini_nus = bio->ini_nus;

  if (setup_exchange_flag)
  {
    setup_exchange(kinetics->grid, kinetics->subgrid.get_box(), { xbcflag == 0, ybcflag == 0, zbcflag == 0 });
    setup_exchange_flag = false;
  }

  DecompGrid<FixKineticsDiffusion>::exchange();

  for (int i = 1; i <= nnus; i++) {
    if (bio->nustate[i] == 0 && !nuConv[i]) {
      if (unit == MOL) {
        xbcm = ini_nus[i][1] * 1000;
        xbcp = ini_nus[i][2] * 1000;
        ybcm = ini_nus[i][3] * 1000;
        ybcp = ini_nus[i][4] * 1000;
        zbcm = ini_nus[i][5] * 1000;
        zbcp = ini_nus[i][6] * 1000;
      } else {
        xbcm = ini_nus[i][1];
        xbcp = ini_nus[i][2];
        ybcm = ini_nus[i][3];
        ybcp = ini_nus[i][4];
        zbcm = ini_nus[i][5];
        zbcp = ini_nus[i][6];
      }
      // copy current concentrations
      for (int grid = 0; grid < snxx_yy_zz; grid++) {
        nuprev[i][grid] = nugrid[i][grid];
      }

      int count = 0;
      // solve diffusion and reaction
      for (int grid = 0; grid < snxx_yy_zz; grid++) {
        // transform nXYZ index to nuR index
        if (ghost[grid] == REGULAR) {
          count++;
          int ind = get_index(grid);
          double nur_ = (unit == KG) ? nur[i][ind] : nur[i][ind] * 1000;
          double diff_coeff;

          if (dcflag) diff_coeff = grid_diff_coeff[i][grid];
          else diff_coeff = bio->diff_coeff[i];

          compute_flux(diff_coeff, nugrid[i][grid], nuprev[i], nur_, grid, ind);

          if (nugrid[i][grid] > 0) {
            (unit == KG) ? (nus[i][ind] = nugrid[i][grid]) : (nus[i][ind] = nugrid[i][grid] / 1000);
          } else {
            nugrid[i][grid] = 1e-20;
            nus[i][ind] = 1e-20;
          }
        }
        if (ghost[grid] == BOUNDARY) {
          double nubs_ = nubs[i];
          
          if(unit == MOL) nubs_ *= 1000;
          compute_bc(nugrid[i][grid], nuprev[i], grid, nubs_);
        }
      }
    }
  }

  int nrequests = 0;
  for (int i = 1; i <= nnus; i++) {
    // checking if is liquid
    if (bio->nustate[i] == 0 && !nuConv[i]) {
      // check convergence criteria
      nuConv[i] = false;
      double max_residual = 0;

      for (int grid = 0; grid < snxx_yy_zz; grid++) {
        if (ghost[grid] == REGULAR) {
          double rate = nugrid[i][grid];
          double prevRate = nuprev[i][grid];
          double residual = fabs((nugrid[i][grid] - nuprev[i][grid]) / nuprev[i][grid]);

          if (residual > max_residual)
            max_residual = residual;
        }
      }

      if (max_residual < tol) nuConv[i] = true;

#if MPI_VERSION >= 3
      MPI_Iallreduce(MPI_IN_PLACE, &nuConv[i], 1, MPI_INT, MPI_BAND, world, &requests[nrequests++]);
#else
      MPI_Allreduce(MPI_IN_PLACE, &nuConv[i], 1, MPI_INT, MPI_BAND, world);
#endif
    }
  }

#if MPI_VERSION >= 3
  MPI_Waitall(nrequests, requests, MPI_STATUS_IGNORE);
#endif

  return nuConv;
}

/* ----------------------------------------------------------------------
 Update grid concentration
  ------------------------------------------------------------------------- */

void FixKineticsDiffusion::update_grids() {
  if (kinetics->blayer < 0)
    return;

  bzhi = kinetics->bnz * stepz;
  if (kinetics->bgrids == 0)
    snxx_yy_zz = 0;
  else
    snxx_yy_zz = snxx * snyy * (MIN(kinetics->subn[2], MAX(0, kinetics->bnz - kinetics->subnlo[2])) + 2);

  for (int grid = 0; grid < snxx_yy_zz; grid++) {
    ghost[grid] = REGULAR;
    if (xgrid[grid][0] < kinetics->sublo[0] || xgrid[grid][1] < kinetics->sublo[1] || xgrid[grid][2] < kinetics->sublo[2]
        || xgrid[grid][0] > kinetics->subhi[0] || xgrid[grid][1] > kinetics->subhi[1] || xgrid[grid][2] > kinetics->subhi[2])
      ghost[grid] = GHOST;
    if (xgrid[grid][0] < xlo || xgrid[grid][1] < ylo || xgrid[grid][2] < zlo
        || xgrid[grid][0] > xhi || xgrid[grid][1] > yhi || xgrid[grid][2] > bzhi)
      ghost[grid] = BOUNDARY;
  }

  if (!bulkflag) return;
  double *nubs = kinetics->nubs;

  for (int grid = 0; grid < snxx * snyy * snzz; grid++) {
    if (xgrid[grid][2] > MIN(bzhi, kinetics->subhi[2])) {
      for (int nu = 1; nu <= bio->nnu; nu++) {
        if (bio->nustate[nu] != 0)
          continue;

        if (kinetics->nubs[nu] == 1e-20) {
          nugrid[nu][grid] = nubs[nu];
        } else {
          nugrid[nu][grid] = (unit == KG) ? nubs[nu] : nubs[nu] * 1000;
        }

        int ind = get_index(grid);
        if (ind != -1) {
          kinetics->nus[nu][ind] = nubs[nu];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
 Update dynamic diffusion coeff
 ------------------------------------------------------------------------- */
void FixKineticsDiffusion::update_diff_coeff() {

  for (int i = 1; i <= bio->nnu; i++) {
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      int ind = get_index(grid);

      if (ind != -1 && kinetics->xdensity[0][ind]) {
        if(dcflag == 1) {
          grid_diff_coeff[i][grid] = bio->diff_coeff[i] * (1 - (0.43 * pow(kinetics->xdensity[0][ind]/vol,0.92)) / (11.19 + 0.27 * pow(kinetics->xdensity[0][ind]/vol,0.99)));
        } else if (dcflag == 2) {
          grid_diff_coeff[i][grid] = bio->diff_coeff[i] * 0.8;
        }
        continue;
      }

      grid_diff_coeff[i][grid] = bio->diff_coeff[i];
    }
  }
}

/* ----------------------------------------------------------------------
 Initialize grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::init_grid() {
  double *nubs = kinetics->nubs;
  double **ini_nus = bio->ini_nus;

  bzhi = kinetics->bnz * stepz;
  if (kinetics->bgrids == 0)
    snxx_yy_zz = 0;
  else
    snxx_yy_zz = snxx * snyy * (MIN(kinetics->subn[2], MAX(0, kinetics->bnz - kinetics->subnlo[2])) + 2);

  double i, j, k;
  int cell = 0;
  for (k = kinetics->sublo[2] - (stepz / 2); k < kinetics->subhi[2] + stepz; k += stepz) {
    for (j = kinetics->sublo[1] - (stepy / 2); j < kinetics->subhi[1] + stepy; j += stepy) {
      for (i = kinetics->sublo[0] - (stepx / 2); i < kinetics->subhi[0] + stepx; i += stepx) {
        xgrid[cell][0] = i;
        xgrid[cell][1] = j;
        xgrid[cell][2] = k;
        //Initialise concentration values for ghost and std grids
        for (int nu = 1; nu <= bio->nnu; nu++) {
          ghost[cell] = REGULAR;
          nugrid[nu][cell] = ini_nus[nu][0];
          if (i < kinetics->sublo[0] || i > kinetics->subhi[0]
              || j < kinetics->sublo[1] || j > kinetics->subhi[1]
              || k < kinetics->sublo[2] || k > kinetics->subhi[2]) {
            ghost[cell] = GHOST;
          } 
          if (i < xlo) {
            ghost[cell] = BOUNDARY;
            nugrid[nu][cell] = ini_nus[nu][1];
          } else if (i > xhi) {
            ghost[cell] = BOUNDARY;
            nugrid[nu][cell] = ini_nus[nu][2];
          } else if (j < ylo) {
            ghost[cell] = BOUNDARY;
            nugrid[nu][cell] = ini_nus[nu][3];
          } else if (j > yhi) {
            ghost[cell] = BOUNDARY;
            nugrid[nu][cell] = ini_nus[nu][4];
          } else if (k < zlo) {
            ghost[cell] = BOUNDARY;
            nugrid[nu][cell] = ini_nus[nu][5];
          } else if (k > MIN(bzhi, zhi)) {
            ghost[cell] = BOUNDARY;
            nugrid[nu][cell] = ini_nus[nu][6];
          }
          if (unit == MOL)
            nugrid[nu][cell] = nugrid[nu][cell] * 1000;
          if (cell == 0)
            nubs[nu] = ini_nus[nu][6];
        }
        cell++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 Mass balance in bulk liquid
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bulk() {
  if (!bulkflag) return;

  double vol = stepx * stepy * stepz;
  double **ini_nus = bio->ini_nus;
  double **nur = kinetics->nur;

  for (int nu = 1; nu <= bio->nnu; nu++) {
    double nubs_ = kinetics->nubs[nu];
    // convert nubs from l to m3
    if (unit == MOL) nubs_ *= 1000;
    // the concentration of o2 in the bulk liquid is kept constant by aeration
    if (bio->nustate[nu] != 0 || !strcmp(bio->nuname[nu], "o2") || !strcmp(bio->nuname[nu], "co2")
        || !strcmp(bio->nuname[nu], "na") || !strcmp(bio->nuname[nu], "cl"))
      continue;

    double sumR = 0;
    double global_sumR = 0;
    // unit in m3
    double inibc = (unit == KG) ? ini_nus[nu][6] : ini_nus[nu][6] * 1000;
    // sum up consumption
    for (int i = 0; i < kinetics->bgrids; i++) {
      (unit == KG) ? (sumR += nur[nu][i]) : (sumR += nur[nu][i] * 1000);
    }

    MPI_Allreduce(&sumR, &global_sumR, 1, MPI_DOUBLE, MPI_SUM, world);

    double dt = update->dt * kinetics->nevery;
    // solve for the mass balance in bulk liquid
    nubs_ = nubs_ + ((q / rvol) * (inibc - nubs_) + ((af * global_sumR * vol) / (rvol * yhi * xhi))) * dt;

    if (unit == MOL) nubs_ /= 1000;
    if (nubs_ < 0) nubs_ = 1e-20;

    kinetics->nubs[nu] = nubs_;
  }
}

/* ----------------------------------------------------------------------
 get index of non-ghost mesh grid
 ------------------------------------------------------------------------- */

int FixKineticsDiffusion::get_index(int grid) {
  int ind;
  int k = grid / snxx_yy;
  int r = grid - k * snxx_yy;
  int j = r / snxx;
  int i = r - j * snxx;
  //if (comm->me == 0) fprintf(screen, "%d -> %d %d %d %d\n", grid, i, j, k, (i - 1) + kinetics->subn[0] * (j - 1) + kinetics->subn[0] * kinetics->subn[1] * (k - 1));
  ind = (i - 1) + kinetics->subn[0] * (j - 1) + kinetics->subn[0] * kinetics->subn[1] * (k - 1);

  if (ind < 0 || ind >= snxyz) return -1;
  else return ind;
}

/* ----------------------------------------------------------------------
 update concentration for ghost grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bc(double &nuCell, double *nuPrev, int grid, double bulk) {
  //for nx = ny = nz = 1 grids
  //18 19 20        21 22 23       24 25 26
  //9  10 11        12 13 14       15 16 17
  //0   1  2         3  4  5        6  7  8
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - snxx;  // y direction
  int fwd = grid + snxx;  // y direction
  int down = grid - snxx * snyy; // z direction
  int up = grid + snxx * snyy;  // z dirction

  // assign values to the ghost-grids according to the boundary conditions.
  // If ghostcells are Neu then take the values equal from the adjacent cells.
  // if ghostcells are dirich then take the values equal to negative of the adjacent cells.
  // if ghostcells are mixed then zlo ghost cells are nuemann, zhi ghost cells are dirichlet, other four surfaces are periodic BC.
  // low-z surface
  if (xgrid[grid][2] < zlo && ghost[up] == REGULAR) {
    //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
    if (zbcflag == PP && kinetics->nz == kinetics->subn[2]) {
      int zhiGrid = grid + snxx * snyy * nz;
      nuCell = nuPrev[zhiGrid];
    } else if (zbcflag == DD) {
      nuCell = 2 * zbcm - nuPrev[up];
    } else if (zbcflag == ND) {
      nuCell = nuPrev[up];
    } else if (zbcflag == NN) {
      nuCell = nuPrev[up];
    } else if (zbcflag == DN) {
      nuCell = 2 * zbcm - nuPrev[up];
    }
  }
  // high-z surface
  else if (xgrid[grid][2] > bzhi && ghost[down] == REGULAR) {
    if (zbcflag == PP && kinetics->nz == kinetics->subn[2]) {
      int zloGrid = grid - snxx * snyy * nz;
      nuCell = nuPrev[zloGrid];
    } else if (zbcflag == DD) {
      nuCell = 2 * bulk - nuPrev[down];
    } else if (zbcflag == ND) {
      nuCell = 2 * bulk - nuPrev[down];
    } else if (zbcflag == NN) {
      nuCell = nuPrev[down];
    } else if (zbcflag == DN) {
      nuCell = nuPrev[down];
    }
  }
  // low-y surface
  else if (xgrid[grid][1] < ylo && ghost[fwd] == REGULAR) {
    if (ybcflag == PP && kinetics->ny == kinetics->subn[1]) {
      int yhiGrid = grid + snxx * ny;
      nuCell = nuPrev[yhiGrid];
    } else if (ybcflag == DD) {
      nuCell = 2 * ybcm - nuPrev[fwd];
    } else if (ybcflag == ND) {
      nuCell = nuPrev[fwd];
    } else if (ybcflag == NN) {
      nuCell = nuPrev[fwd];
    } else if (ybcflag == DN) {
      nuCell = 2 * ybcm - nuPrev[fwd];
    }
  }
  // high-y surface
  else if (xgrid[grid][1] > yhi && ghost[bwd] == REGULAR) {
    if (ybcflag == PP && kinetics->ny == kinetics->subn[1]) {
      int yloGrid = grid - snxx * ny;
      nuCell = nuPrev[yloGrid];
    } else if (ybcflag == DD) {
      nuCell = 2 * ybcp - nuPrev[bwd];
    } else if (ybcflag == ND) {
      nuCell = 2 * ybcp - nuPrev[bwd];
    } else if (ybcflag == NN) {
      nuCell = nuPrev[bwd];
    } else if (ybcflag == DN) {
      nuCell = nuPrev[bwd];
    }
  }
  // low-x surface
  else if (xgrid[grid][0] < xlo && ghost[rhs] == REGULAR) {
    if (xbcflag == PP && kinetics->nx == kinetics->subn[0]) {
      int xhiGrid = grid + nx;
      nuCell = nuPrev[xhiGrid];
    } else if (xbcflag == DD) {
      nuCell = 2 * xbcm - nuPrev[rhs];
    } else if (xbcflag == ND) {
      nuCell = nuPrev[rhs];
    } else if (xbcflag == NN) {
      nuCell = nuPrev[rhs];
    } else if (xbcflag == DN) {
      nuCell = 2 * xbcm - nuPrev[rhs];
    }
  }
  // high-x surface
  else if (xgrid[grid][0] > xhi && ghost[lhs] == REGULAR) {
    if (xbcflag == PP && kinetics->nx == kinetics->subn[0]) {
      int xloGrid = grid - nx;
      nuCell = nuPrev[xloGrid];
    } else if (xbcflag == DD) {
      nuCell = 2 * xbcp - nuPrev[lhs];
    } else if (xbcflag == ND) {
      nuCell = 2 * xbcp - nuPrev[lhs];
    } else if (xbcflag == NN) {
      nuCell = nuPrev[lhs];
    } else if (xbcflag == DN) {
      nuCell = nuPrev[lhs];
    }
  }
}

/* ----------------------------------------------------------------------
 update concentration for non-ghost grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_flux(double cellDNu, double &nuCell, double *nuPrev, double rateNu, int grid, int ind) {
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - snxx;  // y direction
  int fwd = grid + snxx;  // y direction
  int down = grid - snxx * snyy; // z direction
  int up = grid + snxx * snyy;  // z direction

  double jRight = cellDNu * (nuPrev[rhs] - nuPrev[grid]) / stepx;
  double jLeft = cellDNu * (nuPrev[grid] - nuPrev[lhs]) / stepx;
  double jX = (jRight - jLeft) / stepx;

  double jForward = cellDNu * (nuPrev[fwd] - nuPrev[grid]) / stepy;
  double jBackward = cellDNu * (nuPrev[grid] - nuPrev[bwd]) / stepy;
  double jY = (jForward - jBackward) / stepy;

  double jUp = cellDNu * (nuPrev[up] - nuPrev[grid]) / stepz;
  double jDown = cellDNu * (nuPrev[grid] - nuPrev[down]) / stepz;
  double jZ = (jUp - jDown) / stepz;

  double res = 0;
  double shear = 0;

  res = jX + jY + jZ + rateNu;

  // Adding fluxes in all the directions and the uptake rate (RHS side of the equation)
  if (dragflag) {
    double uX = kinetics->fv[0][ind] * (nuPrev[rhs] - nuPrev[lhs]) / (2 * stepx);
    double uY = kinetics->fv[1][ind] * (nuPrev[fwd] - nuPrev[bwd]) / (2 * stepy);
    double uZ = kinetics->fv[2][ind] * (nuPrev[up] - nuPrev[down]) / (2 * stepz);

    res -= (uX + uY + uZ);
    //printf("grid = %i, ind = %i, %e %e %e \n", grid, ind, kinetics->fV[0][ind], kinetics->fV[1][ind], kinetics->fV[2][ind]);
  } else if (shearflag) {
    int hgrid = grid;
    int deep = 0;

    while (ghost[hgrid] == REGULAR) {
      hgrid = hgrid - snxx * snyy;
      deep++;
    }

    shear = srate * (deep * stepz - stepz / 2) * (nuPrev[rhs] - nuPrev[lhs]) / (2 * stepz);

    res -= shear;
  }

  //Updating the value: Ratesub*diffT + nuCell[cell](previous)
  nuCell = nuPrev[grid] + res * diff_dt;
}

/* ----------------------------------------------------------------------
 Compare double values for equality
 ------------------------------------------------------------------------- */

bool FixKineticsDiffusion::is_equal(double a, double b, double c) {
  double epsilon = 1e-10;
  if ((fabs(a - b) > epsilon) || (fabs(a - c) > epsilon) || (fabs(b - c) > epsilon))
    return false;

  return true;
}

/* ----------------------------------------------------------------------
 Manually update reaction if none of the surface is using dirichlet bc
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::update_nus() {
  if ((xbcflag == PP || xbcflag == NN) && (ybcflag == PP || ybcflag == NN) && (zbcflag == PP || zbcflag == NN)) {
    double **nus = kinetics->nus;
    double **nur = kinetics->nur;

    for (int nu = 1; nu < bio->nnu + 1; nu++) {
      if (bio->nustate[nu] != 0)
        continue;

      for (int grid = 0; grid < snxx_yy_zz; grid++) {
        // transform nXYZ index to nuR index
        if (!ghost[grid]) {
          int ind = get_index(grid);
          int dt = update->dt * kinetics->nevery;
          double r = (unit == KG) ? nur[nu][ind] * dt : nur[nu][ind] * dt * 1000;

          nugrid[nu][grid] += r;

          if (nugrid[nu][grid] <= 0)
            nugrid[nu][grid] = 1e-20;

          nus[nu][ind] = (unit == KG) ? nugrid[nu][grid] : nugrid[nu][grid] / 1000;
        }
      }
    }
  }
}

int FixKineticsDiffusion::get_elem_per_cell() const {
  return bio->nnu;
}

void FixKineticsDiffusion::resize(const Subgrid<double, 3> &subgrid) {
  int nnus = bio->nnu;
  snxx_yy_zz = subgrid.cell_count();
  xgrid = memory->grow(xgrid, snxx_yy_zz, 3, "diffusion:xGrid");
  nugrid = memory->grow(nugrid, nnus + 1, snxx_yy_zz, "diffusion:nuGrid");
  nuprev = memory->grow(nuprev, nnus + 1, snxx_yy_zz, "diffusion:nuPrev");
  ghost = memory->grow(ghost, snxx_yy_zz, "diffusion:ghost");
}

void FixKineticsDiffusion::migrate(const Grid<double, 3> &grid, const Box<int, 3> &from, const Box<int, 3> &to) {
  DecompGrid<FixKineticsDiffusion>::migrate(grid, from, to, extend(from), extend(to));
  Subgrid<double, 3> subgrid(grid, to);
  setup_exchange(grid, to, { xbcflag == 0, ybcflag == 0, zbcflag == 0 });
  Subgrid<double, 3> extended(kinetics->grid, extend(to));
  auto cell_centers = extended.get_cell_centers();
  for (int i = 0; i < extended.cell_count(); i++) {
    xgrid[i][0] = cell_centers[i][0];
    xgrid[i][1] = cell_centers[i][1];
    xgrid[i][2] = cell_centers[i][2];
  }
  nx = subgrid.get_dimensions()[0];
  ny = subgrid.get_dimensions()[1];
  nz = subgrid.get_dimensions()[2];
  snxx = extended.get_dimensions()[0];
  snyy = extended.get_dimensions()[1];
  snzz = extended.get_dimensions()[2];
  for (int grid = 0; grid < snxx_yy_zz; grid++) {
    if (xgrid[grid][0] < kinetics->sublo[0] || xgrid[grid][1] < kinetics->sublo[1]
        || xgrid[grid][2] < kinetics->sublo[2] || xgrid[grid][0] > kinetics->subhi[0]
        || xgrid[grid][1] > kinetics->subhi[1] || xgrid[grid][2] > kinetics->subhi[2])
      ghost[grid] = GHOST;
    else
      ghost[grid] = REGULAR;
  }
  double *nubs = kinetics->nubs;
  for (int i = 1; i <= bio->nnu; i++) {
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      if (ghost[grid] == BOUNDARY) {
        compute_bc(nugrid[i][grid], nugrid[i], grid, nubs[i]);
      }
    }
  }
}
