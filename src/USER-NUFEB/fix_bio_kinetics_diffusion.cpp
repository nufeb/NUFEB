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

  if (narg < 4)
    error->all(FLERR, "Not enough arguments in fix diffusion command");

  unit = KG;
  shearflag = dragflag = 0;
  bulkflag = 0;
  srate = 0;
  dcflag = 0;
  closed_system = 0;
  dcratio = 0.8;
  xpbcflag = ypbcflag = zpbcflag = 0;

  xgrid = NULL;
  nugrid = NULL;
  nuprev = NULL;
  grid_type = NULL;
  nuclose = NULL;
  grid_diff_coeff = NULL;
  nupenult = NULL;

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[3 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3 + i][2]);
  }

  int iarg = 4;
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
    } else if (strcmp(arg[iarg], "dcratio") == 0) {
      dcratio = force->numeric(FLERR, arg[iarg + 1]);
      iarg += 2;
      if (dcratio <= 0)
	error->all(FLERR, "Illegal fix kinetics/diffusion command: dcratio");
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
  memory->destroy(grid_type);
  memory->destroy(nuclose);
  if (closed_system) memory->destroy(nupenult);

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
  int energy = 0;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "shear") == 0) {
      shearflag = 1;
    } else if (strcmp(modify->fix[j]->style, "fdrag") == 0) {
      dragflag = 1;
    } else if (strcmp(modify->fix[j]->style, "kinetics/growth/energy") == 0) {
      energy = 1;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for running iBM simulation");

  bio = kinetics->bio;
  tol = input->variable->compute_equal(ivar[0]);

  if (energy) unit = MOL;

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

  snxyz = kinetics->subn[0] * kinetics->subn[1] * kinetics->subn[2];

  snxx = kinetics->subn[0] + 2;
  snyy = kinetics->subn[1] + 2;
  snzz = kinetics->subn[2] + 2;

  snxx_yy = snxx * snyy;
  snxx_yy_zz = snxx_yy * snzz;
  
  int nnu = bio->nnu;

  //inlet concentration and maximum boundary condition conc value
  if (xgrid == NULL) {
    xgrid = memory->create(xgrid, snxx_yy_zz, 3, "diffusion:xgrid");
    nugrid = memory->create(nugrid, nnu+1, snxx_yy_zz, "diffusion:nugrid");
    nuprev = memory->create(nuprev, nnu+1, snxx_yy_zz, "diffusion:nuprev");
    grid_type = memory->create(grid_type, snxx_yy_zz, "diffusion:ghost");
    grid_diff_coeff = memory->create(grid_diff_coeff, nnu+1, snxx_yy_zz, "diffusion:grid_diff_coeff");
    nuclose = memory->create(nuclose, nnu+1, "diffusion:close_nu");
    init_grid();
    init_setting();

    if (closed_system) nupenult = memory->create(nupenult, nnu+1, snxx_yy_zz, "diffusion:nupenult");

    requests = new MPI_Request[MAX(2 * comm->nprocs, nnu + 1)];
  }


  setup_exchange(kinetics->grid, kinetics->subgrid.get_box(), { xpbcflag == 1, ypbcflag == 1, zpbcflag == 1 });

  if (closed_system){
    if (comm->me == 0 && logfile)
      fprintf(logfile, "Model is defined as a closed system. \n");
    if (comm->me == 0 && screen)
      fprintf(screen, "Model is defined as a closed system. \n");
    if (kinetics->blayer >= 0)
      error->all(FLERR, "Cannot define boundary layer in a closed system");
  }
}


/* ----------------------------------------------------------------------
 Initialise diffusion system settings
 ------------------------------------------------------------------------- */
void FixKineticsDiffusion::init_setting() {
  for (int nu = 1; nu <= bio->nnu; nu++) {

    if (bio->nubc[nu][0] == PP)
      xpbcflag = 1;
    if (bio->nubc[nu][1] == PP)
      ypbcflag = 1;
    if (bio->nubc[nu][2] == PP)
      zpbcflag = 1;

    if ((bio->nubc[nu][0] == PP || bio->nubc[nu][0] == NN)
	&& (bio->nubc[nu][1] == PP || bio->nubc[nu][1] == NN)
	&& (bio->nubc[nu][2] == PP || bio->nubc[nu][2] == NN)){
      nuclose[nu] = 1;
      closed_system = 1;
    } else {
      nuclose[nu] = 0;
    }

    if (unit == MOL) {
      init_nusbc = bio->init_nus[nu][1] * M2L;
    } else {
      init_nusbc = bio->init_nus[nu][1];
    }
  }
}

/* ----------------------------------------------------------------------
 Initialise grid attributes
 ------------------------------------------------------------------------- */
void FixKineticsDiffusion::init_grid() {
  double *nubs = kinetics->nubs;
  double **init_nus = bio->init_nus;

  bzhi = kinetics->bnz * stepz;
  if (kinetics->bgrids == 0)
    snxx_yy_zz = 0;
  else
    snxx_yy_zz = snxx * snyy * (MIN(kinetics->subn[2], MAX(0, kinetics->bnz - kinetics->subnlo[2])) + 2);

  double i, j, k;
  int grid = 0;
  for (k = kinetics->sublo[2] - (stepz / 2); k < kinetics->subhi[2] + stepz; k += stepz) {
    for (j = kinetics->sublo[1] - (stepy / 2); j < kinetics->subhi[1] + stepy; j += stepy) {
      for (i = kinetics->sublo[0] - (stepx / 2); i < kinetics->subhi[0] + stepx; i += stepx) {
        xgrid[grid][0] = i;
        xgrid[grid][1] = j;
        xgrid[grid][2] = k;
        //Initialise concentration values for ghost, boundary, and regular grids
        grid_type[grid] = REGULAR;
        if (i < kinetics->sublo[0] || i > kinetics->subhi[0]
            || j < kinetics->sublo[1] || j > kinetics->subhi[1]
            || k < kinetics->sublo[2] || k > kinetics->subhi[2]) {
          grid_type[grid] = GHOST;
        }
        if (i < xlo || i > xhi || j < ylo || j > yhi || k < zlo || k > MIN(bzhi, zhi)) {
          grid_type[grid] = BOUNDARY;
        }

        for (int nu = 1; nu <= bio->nnu; nu++) {
          nugrid[nu][grid] = init_nus[nu][0];

          if (grid_type[grid] == BOUNDARY) {
	    nugrid[nu][grid] = init_nus[nu][1];
          }

          if (unit == MOL)
            nugrid[nu][grid] = nugrid[nu][grid] * M2L;

          if (grid == 0)
            nubs[nu] = init_nus[nu][1];
        }
        grid++;
      }
    }
  }
}


/* ----------------------------------------------------------------------
 Solve diffusion and reaction
 ------------------------------------------------------------------------- */

int *FixKineticsDiffusion::diffusion(int *nuconv, int iter, double diff_dt) {
  int nnu = bio->nnu;
  this->diff_dt = diff_dt;
  double **nur = kinetics->nur;
  double **nus = kinetics->nus;
  double *nubs = kinetics->nubs;
  double **init_nus = bio->init_nus;
  int **nubc = bio->nubc;
  
  if (setup_exchange_flag)
  {
    setup_exchange(kinetics->grid, kinetics->subgrid.get_box(), { xpbcflag == 1, ypbcflag == 1, zpbcflag == 1 });
    setup_exchange_flag = false;
  }

  DecompGrid<FixKineticsDiffusion>::exchange();

  for (int nu = 1; nu <= nnu; nu++) {
    if (bio->nustate[nu] || nuconv[nu]) continue;
      // copy current concentrations
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      if (nuclose[nu] && iter > 1) nupenult[nu][grid] = nuprev[nu][grid];
      nuprev[nu][grid] = nugrid[nu][grid];
    }

    // solve diffusion and reaction
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      // transform nXYZ index to nuR index
       if (grid_type[grid] == REGULAR) {
	int ind = get_index(grid);
	double nur_grid = (unit == KG) ? nur[nu][ind] : nur[nu][ind] * M2L;
	double diff_coeff;

	if (dcflag) diff_coeff = grid_diff_coeff[nu][grid];
	else diff_coeff = bio->diff_coeff[nu];

	compute_flux(diff_coeff, nugrid[nu][grid], nuprev[nu], nur_grid, grid, ind);

	if (nugrid[nu][grid] > 0) {
	  (unit == KG) ? (nus[nu][ind] = nugrid[nu][grid]) : (nus[nu][ind] = nugrid[nu][grid] / M2L);
	} else {
	  nugrid[nu][grid] = MIN_NUS;
	  nus[nu][ind] = MIN_NUS;
	}
      } else if (grid_type[grid] == BOUNDARY) {
	double nubs_ = nubs[nu];

	if(unit == MOL) nubs_ *= M2L;
	compute_bc(nugrid[nu][grid], nuprev[nu], grid, nubs_, nubc[nu]);
      }
    }
  }

  check_converge(nuconv);

  return nuconv;
}


/* ----------------------------------------------------------------------
 check convergence criteria
  ------------------------------------------------------------------------- */

void FixKineticsDiffusion::check_converge(int *nuconv) {
  // check convergence criteria
  int nrequests = 0;
  for (int nu = 1; nu <= bio->nnu; nu++) {
    // checking if is liquid
    if (bio->nustate[nu] || nuconv[nu]) continue;
    nuconv[nu] = false;
    double max_residual = 0;

    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      if (grid_type[grid] == REGULAR) {
	double residual = 0;
	if (!nuclose[nu]) {
	  residual = fabs((nugrid[nu][grid] - nuprev[nu][grid]) / nuprev[nu][grid]);
	} else {
	  double res_n = fabs((nugrid[nu][grid]-nuprev[nu][grid])/nuprev[nu][grid]);
	  double res_n2 = fabs((nuprev[nu][grid]-nupenult[nu][grid])/nupenult[nu][grid]);
	  residual= fabs(res_n-res_n2);
	}

	if (residual > max_residual)
	  max_residual = residual;
      }
    }

    if (max_residual < tol) nuconv[nu] = 1;

#if MPI_VERSION >= 3
      MPI_Iallreduce(MPI_IN_PLACE, &nuconv[nu], 1, MPI_INT, MPI_BAND, world, &requests[nrequests++]);
#else
      MPI_Allreduce(MPI_IN_PLACE, &nuconv[nu], 1, MPI_INT, MPI_BAND, world);
#endif
  }

#if MPI_VERSION >= 3
  MPI_Waitall(nrequests, requests, MPI_STATUS_IGNORE);
#endif
}

/* ----------------------------------------------------------------------
 Average nutrient concentration before solving diffusion in closed system.
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::closed_system_init() {
  double **nus = kinetics->nus;

  for (int nu = 1; nu < bio->nnu + 1; nu++) {
    if (bio->nustate[nu] || !nuclose[nu] || !bio->diff_coeff[nu]) continue;
    double sum = 0;
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      if (grid_type[grid] == REGULAR) {
	sum += nugrid[nu][grid];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
    sum /= (nx * ny * nz);
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      nugrid[nu][grid] = sum;
      if (unit == MOL) nugrid[nu][grid] /= M2L;
    }
  }
}

/* ----------------------------------------------------------------------
 Update nutrient concentrations in closed system when closed_flag = 0
 The update triggered after diffusion process, and based on residual value and
 biological timestep
 ------------------------------------------------------------------------- */
void FixKineticsDiffusion::closed_system_scaleup(double dt, double diff_dt) {
  double **nus = kinetics->nus;
  for (int nu = 1; nu <= bio->nnu; nu++) {
    if (bio->nustate[nu] || !nuclose[nu] || !bio->diff_coeff[nu]) continue;

    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      double residual = nugrid[nu][grid] - nuprev[nu][grid];
      nugrid[nu][grid] += residual/diff_dt*dt;
      if (unit == MOL) nugrid[nu][grid] /= M2L;
      if (nugrid[nu][grid] <= 0) {
	nugrid[nu][grid] = MIN_NUS;
      }

      if (grid_type[grid] == REGULAR) {
	int ind = get_index(grid);
	if (ind != -1) nus[nu][ind] = nugrid[nu][grid];
      }
    }
  }
}


/* ----------------------------------------------------------------------
 Update grid attributes
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
    grid_type[grid] = REGULAR;
    if (xgrid[grid][0] < kinetics->sublo[0] || xgrid[grid][1] < kinetics->sublo[1] || xgrid[grid][2] < kinetics->sublo[2]
        || xgrid[grid][0] > kinetics->subhi[0] || xgrid[grid][1] > kinetics->subhi[1] || xgrid[grid][2] > kinetics->subhi[2])
      grid_type[grid] = GHOST;
    if (xgrid[grid][0] < xlo || xgrid[grid][1] < ylo || xgrid[grid][2] < zlo
        || xgrid[grid][0] > xhi || xgrid[grid][1] > yhi || xgrid[grid][2] > bzhi)
      grid_type[grid] = BOUNDARY;
  }

  if (!bulkflag) return;
  double *nubs = kinetics->nubs;

  for (int grid = 0; grid < snxx*snyy*snzz; grid++) {
    if (xgrid[grid][2] > MIN(bzhi, kinetics->subhi[2])) {
      for (int nu = 1; nu <= bio->nnu; nu++) {
        if (bio->nustate[nu] != 0)
          continue;

        if (kinetics->nubs[nu] == MIN_NUS) {
          nugrid[nu][grid] = nubs[nu];
        } else {
          nugrid[nu][grid] = (unit == KG) ? nubs[nu] : nubs[nu] * M2L;
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
          grid_diff_coeff[i][grid] = bio->diff_coeff[i] * dcratio;
        }
        continue;
      }

      grid_diff_coeff[i][grid] = bio->diff_coeff[i];
    }
  }
}


/* ----------------------------------------------------------------------
 Mass balance in bulk liquid
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bulk() {
  if (!bulkflag) return;

  double vol = stepx * stepy * stepz;
  double **init_nus = bio->init_nus;
  double **nur = kinetics->nur;

  for (int nu = 1; nu <= bio->nnu; nu++) {
    double nubs_ = kinetics->nubs[nu];
    // convert nubs from l to m3
    if (unit == MOL) nubs_ *= M2L;
    // the concentration of o2 in the bulk liquid is kept constant by aeration
    if (bio->nustate[nu] != 0 || !strcmp(bio->nuname[nu], "o2") || !strcmp(bio->nuname[nu], "co2")
        || !strcmp(bio->nuname[nu], "na") || !strcmp(bio->nuname[nu], "cl"))
      continue;

    double sumR = 0;
    double global_sumR = 0;
    // unit in m3
    double inibc = (unit == KG) ? init_nus[nu][1] : init_nus[nu][1] * M2L;
    // sum up consumption
    for (int i = 0; i < kinetics->bgrids; i++) {
      (unit == KG) ? (sumR += nur[nu][i]) : (sumR += nur[nu][i] * M2L);
    }

    MPI_Allreduce(&sumR, &global_sumR, 1, MPI_DOUBLE, MPI_SUM, world);

    double dt = update->dt * kinetics->nevery;
    // solve for the mass balance in bulk liquid
    nubs_ = nubs_ + ((q / rvol) * (inibc - nubs_) + ((af * global_sumR * vol) / (rvol * yhi * xhi))) * dt;

    if (unit == MOL) nubs_ /= M2L;
    if (nubs_ < 0) nubs_ = MIN_NUS;

    kinetics->nubs[nu] = nubs_;
  }
}

/* ----------------------------------------------------------------------
 get index of non-ghost mesh grid
 ------------------------------------------------------------------------- */

int FixKineticsDiffusion::get_index(int grid) {
  int ind;
  int k = grid/snxx_yy;
  int r = grid-k*snxx_yy;
  int j = r/snxx;
  int i = r-j*snxx;
  //if (comm->me == 0) fprintf(screen, "%d -> %d %d %d %d\n", grid, i, j, k, (i - 1) + kinetics->subn[0] * (j - 1) + kinetics->subn[0] * kinetics->subn[1] * (k - 1));
  ind = (i-1) + kinetics->subn[0] * (j-1) + kinetics->subn[0] * kinetics->subn[1] * (k-1);

  if (ind < 0 || ind >= snxyz) return -1;
  else return ind;
}

/* ----------------------------------------------------------------------
 update concentration for ghost grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bc(double &nus_grid, double *nus_prev, int grid, double bulk, int *nubc) {
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
  // low-z surface
  if (xgrid[grid][2] < zlo && grid_type[up] == REGULAR) {
    //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
    if (nubc[2] == PP && kinetics->nz == kinetics->subn[2]) {
      int zhiGrid = grid + snxx * snyy * nz;
      nus_grid = nus_prev[zhiGrid];
    } else if (nubc[2] == DD) {
      nus_grid = 2 * init_nusbc - nus_prev[up];
    } else if (nubc[2] == ND) {
      nus_grid = nus_prev[up];
    } else if (nubc[2] == NN) {
      nus_grid = nus_prev[up];
    } else if (nubc[2] == DN) {
      nus_grid = 2 * init_nusbc - nus_prev[up];
    }
  }
  // high-z surface
  else if (xgrid[grid][2] > bzhi && grid_type[down] == REGULAR) {
    if (nubc[2] == PP && kinetics->nz == kinetics->subn[2]) {
      int zloGrid = grid - snxx * snyy * nz;
      nus_grid = nus_prev[zloGrid];
    } else if (nubc[2] == DD) {
      nus_grid = 2 * bulk - nus_prev[down];
    } else if (nubc[2] == ND) {
      nus_grid = 2 * bulk - nus_prev[down];
    } else if (nubc[2] == NN) {
      nus_grid = nus_prev[down];
    } else if (nubc[2] == DN) {
      nus_grid = nus_prev[down];
    }
  }
  // low-y surface
  else if (xgrid[grid][1] < ylo && grid_type[fwd] == REGULAR) {
    if (nubc[1] == PP && kinetics->ny == kinetics->subn[1]) {
      int yhiGrid = grid + snxx * ny;
      nus_grid = nus_prev[yhiGrid];
    } else if (nubc[1] == DD) {
      nus_grid = 2 * init_nusbc - nus_prev[fwd];
    } else if (nubc[1] == ND) {
      nus_grid = nus_prev[fwd];
    } else if (nubc[1] == NN) {
      nus_grid = nus_prev[fwd];
    } else if (nubc[1] == DN) {
      nus_grid = 2 * init_nusbc - nus_prev[fwd];
    }
  }
  // high-y surface
  else if (xgrid[grid][1] > yhi && grid_type[bwd] == REGULAR) {
    if (nubc[1] == PP && kinetics->ny == kinetics->subn[1]) {
      int yloGrid = grid - snxx * ny;
      nus_grid = nus_prev[yloGrid];
    } else if (nubc[1] == DD) {
      nus_grid = 2 * init_nusbc - nus_prev[bwd];
    } else if (nubc[1] == ND) {
      nus_grid = 2 * init_nusbc - nus_prev[bwd];
    } else if (nubc[1] == NN) {
      nus_grid = nus_prev[bwd];
    } else if (nubc[1] == DN) {
      nus_grid = nus_prev[bwd];
    }
  }
  // low-x surface
  else if (xgrid[grid][0] < xlo && grid_type[rhs] == REGULAR) {
    if (nubc[0] == PP && kinetics->nx == kinetics->subn[0]) {
      int xhiGrid = grid + nx;
      nus_grid = nus_prev[xhiGrid];
    } else if (nubc[0] == DD) {
      nus_grid = 2 * init_nusbc - nus_prev[rhs];
    } else if (nubc[0] == ND) {
      nus_grid = nus_prev[rhs];
    } else if (nubc[0] == NN) {
      nus_grid = nus_prev[rhs];
    } else if (nubc[0] == DN) {
      nus_grid = 2 * init_nusbc - nus_prev[rhs];
    }
  }
  // high-x surface
  else if (xgrid[grid][0] > xhi && grid_type[lhs] == REGULAR) {
    if (nubc[0] == PP && kinetics->nx == kinetics->subn[0]) {
      int xloGrid = grid - nx;
      nus_grid = nus_prev[xloGrid];
    } else if (nubc[0] == DD) {
      nus_grid = 2 * init_nusbc - nus_prev[lhs];
    } else if (nubc[0] == ND) {
      nus_grid = 2 * init_nusbc - nus_prev[lhs];
    } else if (nubc[0] == NN) {
      nus_grid = nus_prev[lhs];
    } else if (nubc[0] == DN) {
      nus_grid = nus_prev[lhs];
    }
  }
}

/* ----------------------------------------------------------------------
 update concentration for non-ghost grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_flux(double grid_dc, double &nus_grid, double *nus_prev, double nur_grid, int grid, int ind) {
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - snxx;  // y direction
  int fwd = grid + snxx;  // y direction
  int down = grid - snxx * snyy; // z direction
  int up = grid + snxx * snyy;  // z direction

  double jRight = grid_dc * (nus_prev[rhs] - nus_prev[grid]) / stepx;
  double jLeft = grid_dc * (nus_prev[grid] - nus_prev[lhs]) / stepx;
  double jX = (jRight - jLeft) / stepx;

  double jForward = grid_dc * (nus_prev[fwd] - nus_prev[grid]) / stepy;
  double jBackward = grid_dc * (nus_prev[grid] - nus_prev[bwd]) / stepy;
  double jY = (jForward - jBackward) / stepy;

  double jUp = grid_dc * (nus_prev[up] - nus_prev[grid]) / stepz;
  double jDown = grid_dc * (nus_prev[grid] - nus_prev[down]) / stepz;
  double jZ = (jUp - jDown) / stepz;

  double res = 0;
  double shear = 0;

  res = jX + jY + jZ + nur_grid;

  // Adding fluxes in all the directions and the uptake rate (RHS side of the equation)
  if (dragflag) {
    double uX = kinetics->fv[0][ind] * (nus_prev[rhs] - nus_prev[lhs]) / (2 * stepx);
    double uY = kinetics->fv[1][ind] * (nus_prev[fwd] - nus_prev[bwd]) / (2 * stepy);
    double uZ = kinetics->fv[2][ind] * (nus_prev[up] - nus_prev[down]) / (2 * stepz);

    res -= (uX + uY + uZ);
    //printf("grid = %i, ind = %i, %e %e %e \n", grid, ind, kinetics->fV[0][ind], kinetics->fV[1][ind], kinetics->fV[2][ind]);
  } else if (shearflag) {
    int hgrid = grid;
    int deep = 0;

    while (grid_type[hgrid] == REGULAR) {
      hgrid = hgrid - snxx * snyy;
      deep++;
    }

    shear = srate * (deep * stepz - stepz / 2) * (nus_prev[rhs] - nus_prev[lhs]) / (2 * stepz);

    res -= shear;
  }

  //Updating the value: Ratesub*diffT + nuCell[cell](previous)
  nus_grid = nus_prev[grid] + res * diff_dt;
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
  grid_type = memory->grow(grid_type, snxx_yy_zz, "diffusion:ghost");
}

void FixKineticsDiffusion::migrate(const Grid<double, 3> &grid, const Box<int, 3> &from, const Box<int, 3> &to) {
  DecompGrid<FixKineticsDiffusion>::migrate(grid, from, to, extend(from), extend(to));
  Subgrid<double, 3> subgrid(grid, to);
  setup_exchange(grid, to, { xpbcflag == 1, ypbcflag == 1, zpbcflag == 1 });
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
      grid_type[grid] = GHOST;
    else
      grid_type[grid] = REGULAR;
  }
  double *nubs = kinetics->nubs;
  for (int i = 1; i <= bio->nnu; i++) {
    for (int grid = 0; grid < snxx_yy_zz; grid++) {
      if (grid_type[grid] == BOUNDARY) {
        compute_bc(nugrid[i][grid], nugrid[i], grid, nubs[i], bio->nubc[i]);
      }
    }
  }
}
