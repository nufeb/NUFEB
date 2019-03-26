/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include "fix_bio_kinetics.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm> 

#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"

#include "bio.h"
#include "atom_vec_bio.h"
#include "fix_bio_kinetics_ph.h"
#include "fix_bio_kinetics_thermo.h"
#include "pointers.h"
#include "variable.h"
#include "modify.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "domain.h"
#include "fix_bio_kinetics_diffusion.h"
#include "fix_bio_kinetics_energy.h"
#include "fix_bio_kinetics_monod.h"
#include "comm.h"
#include "fix_bio_fluid.h"
#include "compute.h"
#include "compute_bio_height.h"
#include "compute_bio_rough.h"

#ifdef OUTPUT_GRID
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

#define BUFMIN 1000

/* ---------------------------------------------------------------------- */

FixKinetics::FixKinetics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  bio = avec->bio;

  if (narg < 9)
    error->all(FLERR, "Not enough arguments in fix kinetics command");

  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix kinetics command: calling steps should be positive integer");

  nx = force->inumeric(FLERR, arg[4]);
  ny = force->inumeric(FLERR, arg[5]);
  nz = force->inumeric(FLERR, arg[6]);

  var = new char*[2];
  ivar = new int[2];

  for (int i = 0; i < 2; i++) {
    int n = strlen(&arg[7 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[7 + i][2]);
  }

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

  //set default values
  temp = 298.15;
  rth = 0.0083144;
  demflag = 0;
  niter = -1;
  devery = 1;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temp") == 0) {
      temp = force->numeric(FLERR, arg[iarg + 1]);
      if (temp < 0.0)
        error->all(FLERR, "Illegal fix kinetics command: temp");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rth") == 0) {
      rth = force->numeric(FLERR, arg[iarg + 1]);
      if (rth < 0.0)
        error->all(FLERR, "Illegal fix kinetics command: rth");
      iarg += 2;
    } else if (strcmp(arg[iarg], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (demflag != 1 && demflag != 0)
        error->all(FLERR, "Illegal fix kinetics command: demflag");
      iarg += 2;
    } else if (strcmp(arg[iarg], "niter") == 0) {
      niter = force->inumeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "devery") == 0) {
      devery = force->inumeric(FLERR, arg[iarg + 1]);
      if (devery < 1)
        error->all(FLERR, "Illegal fix kinetics command: devery");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix kinetics command");
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  grid = Grid<double, 3>(Box<double, 3>(domain->boxlo, domain->boxhi), { nx, ny, nz });
  double tmpsublo[3], tmpsubhi[3];
  const double small = 1e-12;
  for (int i = 0; i < 3; i++) {
    tmpsublo[i] = domain->sublo[i] + small;
    tmpsubhi[i] = domain->subhi[i] + small;
  }
  subgrid = Subgrid<double, 3>(grid, Box<double, 3>(tmpsublo, tmpsubhi));

  for (int i = 0; i < 3; i++) {
    // considering that the grid will always have a cubic cell (i.e. stepx = stepy = stepz)
    subnlo[i] = tmpsublo[i] / stepz;
    subnhi[i] = tmpsubhi[i] / stepz;
    sublo[i] = subnlo[i] * stepz;
    subhi[i] = subnhi[i] * stepz;
    subn[i] = subnhi[i] - subnlo[i];
  }

  bnz = subgrid.get_box().upper[2];
  maxheight = domain->boxhi[2];
}

/* ---------------------------------------------------------------------- */

FixKinetics::~FixKinetics() {
  int i;
  for (i = 0; i < 2; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(grid_yield);
  memory->destroy(activity);
  memory->destroy(nur);
  memory->destroy(nus);
  memory->destroy(nubs);
  memory->destroy(gibbs_cata);
  memory->destroy(gibbs_anab);;
  memory->destroy(sh);
  memory->destroy(fv);
  memory->destroy(xdensity);

  delete[] nuconv;
}

/* ---------------------------------------------------------------------- */

int FixKinetics::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init() {
  for (int n = 0; n < 2; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix kinetics is invalid style");
  }

  diff_dt = input->variable->compute_equal(ivar[0]);
  blayer = input->variable->compute_equal(ivar[1]);

  // register fix kinetics with this class
  diffusion = NULL;
  energy = NULL;
  ph = NULL;
  thermo = NULL;
  monod = NULL;
  nufebfoam = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics/growth/energy") == 0) {
      energy = static_cast<FixKineticsEnergy *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/diffusion") == 0) {
      diffusion = static_cast<FixKineticsDiffusion *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/ph") == 0) {
      ph = static_cast<FixKineticsPH *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/thermo") == 0) {
      thermo = static_cast<FixKineticsThermo *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/growth/monod") == 0) {
      monod = static_cast<FixKineticsMonod *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "nufebFoam") == 0) {
      nufebfoam = static_cast<FixFluid *>(lmp->modify->fix[j]);
    }
  }

  if (bio->nnu == 0)
    error->all(FLERR, "fix_kinetics requires # of Nutrients inputs");
  else if (bio->nugibbs_coeff == NULL && energy != NULL)
    error->all(FLERR, "fix_kinetics requires Nutrient Energy inputs");
  else if (bio->ini_nus == NULL)
    error->all(FLERR, "fix_kinetics requires Nutrients inputs");
  if (energy != NULL || thermo != NULL || ph != NULL){
    if (thermo == NULL)
      error->all(FLERR, "fix_kinetics requires fix_kinetics/thermo");
    if (ph == NULL)
      error->all(FLERR, "fix_kinetics requires fix_kinetics/ph");
    if (monod != NULL)
      error->all(FLERR, "kinetics/growth/monod and kinetics/growth/energy cannot be defined at the same time");
  }

  ngrids = subn[0] * subn[1] * subn[2];

  int ntypes = atom->ntypes;
  int nnus = bio->nnu;

  nuconv = new int[nnus + 1]();
  nus = memory->create(nus, nnus + 1, ngrids, "kinetics:nus");
  nur = memory->create(nur, nnus + 1, ngrids, "kinetics:nur");
  nubs = memory->create(nubs, nnus + 1, "kinetics:nubs");
  grid_yield = memory->create(grid_yield, ntypes + 1, ngrids, "kinetic:grid_yield");
  activity = memory->create(activity, nnus + 1, 5, ngrids, "kinetics:activity");
  gibbs_cata = memory->create(gibbs_cata, ntypes + 1, ngrids, "kinetics:gibbs_cata");
  gibbs_anab = memory->create(gibbs_anab, ntypes + 1, ngrids, "kinetics:gibbs_anab");
  sh = memory->create(sh, ngrids, "kinetics:sh");
  fv = memory->create(fv, 3, ngrids, "kinetcis:fv");
  xdensity = memory->create(xdensity, ntypes + 1, ngrids, "kinetics:xdensity");

  // Fitting initial domain decomposition to the grid 
  for (int i = 0; i < comm->procgrid[0]; i++) {
    int n = nx * i * 1.0 / comm->procgrid[0];
    comm->xsplit[i] = (double) n / nx;
  }
  for (int i = 0; i < comm->procgrid[1]; i++) {
    int n = ny * i * 1.0 / comm->procgrid[1];
    comm->ysplit[i] = (double) n / ny;
  }
  for (int i = 0; i < comm->procgrid[2]; i++) {
    int n = nz * i * 1.0 / comm->procgrid[2];
    comm->zsplit[i] = (double) n / nz;
  }
  domain->set_local_box();

  init_param();
  reset_isconv();
  update_bgrids();
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init_param() {
  for (int j = 0; j < ngrids; j++) {
    fv[0][j] = 0;
    fv[1][j] = 0;
    fv[2][j] = 0;

    for (int i = 0; i <= atom->ntypes; i++) {
      grid_yield[i][j] = bio->yield[i];
      gibbs_cata[i][j] = 0;
      gibbs_anab[i][j] = 0;
      xdensity[i][j] = 0;
    }

    for (int i = 0; i <= bio->nnu; i++) {
      nus[i][j] = bio->ini_nus[i][0];
      nur[i][j] = 0;

      activity[i][0][j] = 0;
      activity[i][1][j] = 0;
      activity[i][2][j] = 0;
      activity[i][3][j] = 0;
      activity[i][4][j] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixKinetics::pre_force(int vflag) {
  bool flag = true;

  if (nevery == 0)
    flag = false;
  if (update->ntimestep % nevery)
    flag = false;
  if (nufebfoam != NULL && nufebfoam->demflag)
    flag = false;
  if (demflag)
    flag = false;

  if (flag)
    integration();
}

/* ----------------------------------------------------------------------
 integration loop
 ------------------------------------------------------------------------- */
void FixKinetics::integration() {
  int iteration = 0;
  bool converge = false;
  int nnus = bio->nnu;

  grow_flag = 0;
  update_bgrids();
  update_xdensity();

  // update grid biomass to calculate diffusion coeff
  if (diffusion != NULL && diffusion->dcflag) {
    diffusion->update_diff_coeff();
  }

  while (!converge) {
    converge = true;

    // solve for reaction term, no growth happens here
    if (iteration % devery == 0) {
      reset_nur();
      if (energy != NULL) {
        ph->solve_ph();
        thermo->thermo(diff_dt * devery);
        energy->growth(diff_dt * devery, grow_flag);
      } else if (monod != NULL) {
        monod->growth(diff_dt * devery, grow_flag);
      }
    }

    iteration++;

    // solve for diffusion and advection
    if (diffusion != NULL) {
      nuconv = diffusion->diffusion(nuconv, iteration, diff_dt);
    } else {
      break;
    }

    // check for convergence
    for (int i = 1; i <= nnus; i++) {
      if (!nuconv[i]) {
        converge = false;
        reset_isconv();
        break;
      }
    }

    if (niter > 0 && iteration >= niter)
      converge = true;
  }

  if (comm->me == 0 && logfile)
    fprintf(logfile, "number of iterations: %i \n", iteration);
  if (comm->me == 0 && screen)
    fprintf(screen, "number of iterations: %i \n", iteration);

  grow_flag = 1;
  reset_isconv();
  reset_nur();

  // microbe growth
  if (energy != NULL)
    energy->growth(update->dt * nevery, grow_flag);
  if (monod != NULL)
    monod->growth(update->dt * nevery, grow_flag);

  if (ph != NULL && ph->buffer_flag)
    ph->buffer_ph();

  if (diffusion != NULL) {
    // manually update reaction if none of the surface is using dirichlet BC
    diffusion->update_nus();
    // update concentration in bulk liquid
    diffusion->compute_bulk();
    // update grids
    diffusion->update_grids();
  }

  if (thermo != NULL)
    thermo->thermo(update->dt * nevery);
}

/* ----------------------------------------------------------------------
 get maximum biofilm height
 ------------------------------------------------------------------------- */
double FixKinetics::get_max_height() {
  const int nlocal = atom->nlocal;
  double * const * const x = atom->x;
  double * const r = atom->radius;
  double maxh = 0;

  for (int i = 0; i < nlocal; i++) {
    if ((x[i][2] + r[i]) > maxh)
      maxh = x[i][2] + r[i];
  }

  double global_max;
  MPI_Allreduce(&maxh, &global_max, 1, MPI_DOUBLE, MPI_MAX, world);

  return global_max;
}

/* ----------------------------------------------------------------------
 update grid size based on boundary layer height
 ------------------------------------------------------------------------- */
void FixKinetics::update_bgrids() {
  if (blayer >= 0) {
    maxheight = get_max_height();
    int new_bnz = (int) ((blayer + maxheight) / stepz) + 1;
    if (new_bnz != bnz)
    {
      bnz = new_bnz;
      
      bgrids = subn[0] * subn[1] * MIN(subn[2], MAX(0, bnz - subnlo[2]));
      
      double tmpsublo[3], tmpsubhi[3];
      const double small = 1e-12;
      for (int i = 0; i < 3; i++) {
        tmpsublo[i] = domain->sublo[i] + small;
        if (i == 2)
          tmpsubhi[i] = bnz * stepz;
        else
          tmpsubhi[i] = domain->subhi[i] + small;
      }
      subgrid = Subgrid<double, 3>(grid, Box<double, 3>(tmpsublo, tmpsubhi));
      diffusion->setup_exchange_flag = true;
    }
  } else {
    bgrids = subn[0] * subn[1] * subn[2];
  }
}

/* ----------------------------------------------------------------------
 update biomass density
 ------------------------------------------------------------------------- */
void FixKinetics::update_xdensity() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double vol = stepx * stepy * stepz;

  for (int i = 0; i <= atom->ntypes; i++) {
    for (int j = 0; j < bgrids; j++) {
      xdensity[i][j] = 0;
    }
  }

  for (int i = 0; i < nlocal; i++) {
    int pos = position(i);
    int t = atom->type[i];
    double xmass = atom->rmass[i] / vol;
    xdensity[t][pos] += xmass;
    xdensity[0][pos] += xmass;
  }
}

bool FixKinetics::is_inside(int i) {
  if (atom->x[i][0] < sublo[0] || atom->x[i][0] >= subhi[0] || atom->x[i][1] < sublo[1] || atom->x[i][1] >= subhi[1]
      || atom->x[i][2] < sublo[2] || atom->x[i][2] >= subhi[2])
    return false;
  return true;
}

/* ----------------------------------------------------------------------
 get grid index of atom i
 ------------------------------------------------------------------------- */
int FixKinetics::position(int i) {
  // get index of grid containing i
  int xpos = (atom->x[i][0] - sublo[0]) / stepz;
  int ypos = (atom->x[i][1] - sublo[1]) / stepz;
  int zpos = (atom->x[i][2] - sublo[2]) / stepz;
  int pos = xpos + ypos * subn[0] + zpos * subn[0] * subn[1];

  if (pos >= bgrids) {
    printf("Too big! pos=%d   size = %i\n", pos, bgrids);
  }

  return pos;
}

/* ----------------------------------------------------------------------
 reset nutrient reaction array
 ------------------------------------------------------------------------- */
void FixKinetics::reset_nur() {
  for (int nu = 1; nu < bio->nnu + 1; nu++) {
    for (int j = 0; j < bgrids; j++) {
      nur[nu][j] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
 reset convergence status
 ------------------------------------------------------------------------- */
void FixKinetics::reset_isconv() {
  for (int i = 1; i <= bio->nnu; i++) {
    if (strcmp(bio->nuname[i], "h2o") == 0)
      nuconv[i] = true;
    else if (strcmp(bio->nuname[i], "h") == 0)
      nuconv[i] = true;
    else if (bio->diff_coeff[i] == 0)
      nuconv[i] = true;
    else if (bio->nustate[i] == 1)
      nuconv[i] = true;
    else
      nuconv[i] = false;
  }
}

/* ---------------------------------------------------------------------- */
int FixKinetics::modify_param(int narg, char **arg) {
  if (strcmp(arg[0], "demflag") == 0) {
    if (narg != 2)
      error->all(FLERR, "Illegal fix_modify command");
    demflag = force->inumeric(FLERR, arg[1]);
    if (demflag != 1 && demflag != 0)
      error->all(FLERR, "Illegal fix_modify command: demflag");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixKinetics::get_elem_per_cell() const {
  int result = 2 * bio->nnu; // nuS + nuR
  if (energy) {
    result += 5 * bio->nnu; // qGas + activity
    result += 3 * atom->ntypes; // gYield + DRGCat + DRGAn 
  }
  if (nufebfoam) {
    result += 3; // fV
  }
  return result;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::migrate() {
  Subgrid<double, 3> new_subgrid = Subgrid<double, 3>(grid, Box<double, 3>(domain->sublo, domain->subhi),
      [](double value) {return std::round(value);});
  DecompGrid<FixKinetics>::migrate(grid, subgrid.get_box(), new_subgrid.get_box());
  for (int i = 0; i < 3; i++) {
    subnlo[i] = new_subgrid.get_origin()[i];
    subnhi[i] = subnlo[i] + new_subgrid.get_dimensions()[i];
    sublo[i] = subnlo[i] * stepz;
    subhi[i] = subnhi[i] * stepz;
    subn[i] = subnhi[i] - subnlo[i];
  }
  diffusion->migrate(grid, subgrid.get_box(), new_subgrid.get_box());
  subgrid = new_subgrid;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::resize(const Subgrid<double, 3> &subgrid) {
  int nnus = bio->nnu;
  int ntypes = atom->ntypes;
  ngrids = subgrid.cell_count();
  update_bgrids();
  nus = memory->grow(nus, nnus + 1, ngrids, "kinetics:nus");
  nur = memory->grow(nur, nnus + 1, ngrids, "kinetics:nur");
  if (energy) {
    grid_yield = memory->grow(grid_yield, ntypes + 1, ngrids, "kinetic:grid_yield");
    activity = memory->grow(activity, nnus + 1, 5, ngrids, "kinetics:activity");
    gibbs_cata = memory->grow(gibbs_cata, ntypes + 1, ngrids, "kinetics:gibbs_cata");
    gibbs_anab = memory->grow(gibbs_anab, ntypes + 1, ngrids, "kinetics:gibbs_anab");
    sh = memory->grow(sh, ngrids, "kinetics:sh");
  }
  if (nufebfoam) {
    fv = memory->grow(fv, 3, ngrids, "kinetcis:fV");
  }
  if (monod != NULL)
    monod->grow_subgrid(ngrids);
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->style == "ave_height")
      static_cast<ComputeNufebHeight *>(modify->compute[i])->grow_subgrid();
    else if (modify->compute[i]->style == "roughness")
      static_cast<ComputeNufebRough *>(modify->compute[i])->grow_subgrid();
  }
}
