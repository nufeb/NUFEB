/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_verify.h"

#include <cstdio>
#include <cstring>

#include "atom.h"
#include "atom_vec_bio.h"
#include "bio.h"
#include "error.h"
#include "fix_bio_kinetics.h"
#include "compute_bio_height.h"
#include "fix_bio_kinetics_diffusion.h"
#include "fix_bio_kinetics_monod.h"
#include "force.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include <stdio.h>
#include <math.h>
#include "comm.h"

#include <vector>


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixVerify::FixVerify(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
 // if (narg > 7) error->all(FLERR,"Illegal fix verify command");

  nevery = force->inumeric(FLERR,arg[3]);

  int nkeywords = narg - 4;
  bm1flag = bm2flag = bm3flag = mflag = 0;
  demflag = 0;

  for (int i = 0; i < nkeywords; i++) {
    if (strcmp(arg[4+i], "bm1") == 0) bm1flag = 1;
    else if (strcmp(arg[4+i], "bm2") == 0) bm2flag = 1;
    else if (strcmp(arg[4+i], "mb") == 0) mflag = 1;
    else if (strcmp(arg[4+i], "bm3") == 0) bm3flag = 1;
    else if (strcmp(arg[4+i], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[5+i]);
      if (demflag != 0 && demflag != 1)
        error->all(FLERR, "Illegal fix divide command: demflag");
    }
  }
}

FixVerify::~FixVerify()
{
}

void FixVerify::init()
{
  // get overall concentration
  global_no2 = 0;
  global_nh3 = 0;
  global_pre_no2 = 0;
  global_pre_nh3 = 0;
  global_smass = 0;
  global_pre_smass = 0;

  kinetics = NULL;
  diffusion = NULL;

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/diffusion") == 0) {
      diffusion = static_cast<FixKineticsDiffusion *>(lmp->modify->fix[j]);
    }
  }

  int ncompute = modify->ncompute;

  for (int j = 0; j < ncompute; j++) {
    if (strcmp(modify->compute[j]->style, "ave_height") == 0) {
      cheight = static_cast<ComputeNufebHeight *>(lmp->modify->compute[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required");

  bio = kinetics->bio;
  vol = kinetics->stepx * kinetics->stepy * kinetics->stepz;
 // kinetics->diffusion->bulkflag = 0;

  if (bm2flag || bm1flag || bm3flag) {
    free_particle_list();
  }

  avec = (AtomVecBio *) atom->style_match("bio");
}

/* ---------------------------------------------------------------------- */

int FixVerify::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


void FixVerify::end_of_step() {
  if (update->ntimestep % nevery) return;
  if (demflag) return;

  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;

  catCoeff = bio->cata_coeff;
  anabCoeff = bio->anab_coeff;

  gYield = kinetics->grid_yield;
  nuS = kinetics->nus;
  nnus = bio->nnu;

  if (mflag == 1) nitrogen_mass_balance();
  if (bm1flag == 1) benchmark_one();
  if (bm2flag == 1) benchmark_two();
  if (bm3flag == 1) benchmark_three();
}

/* ----------------------------------------------------------------------
 mass balance check for Nitrogen
 ------------------------------------------------------------------------- */

void FixVerify::nitrogen_mass_balance() {
  double smass, pre_smass;
  double nh3_nitrogen, pre_nh3_nitrogen;
  double no2_nitrogen, pre_no2_nitrogen;
  double no3_nitrogen, pre_no3_nitrogen;
  int ngrids;

  nh3_nitrogen = pre_nh3_nitrogen = 0;
  no2_nitrogen = pre_no2_nitrogen = 0;
  no3_nitrogen = pre_no3_nitrogen = 0;

  smass = 0;

  // get biomass concentration (mol/L)
  for (int i = 0; i < nlocal; i++) {
    int pos = kinetics->position(i);
    double rmassCellVol = atom->rmass[i] / vol;
    rmassCellVol /= 24.6;

    if (pos != -1) smass += rmassCellVol;
  }

  // get overall biamass concentration
  MPI_Allreduce(&smass,&global_smass,1,MPI_DOUBLE,MPI_SUM,world);

  // get nitrogen concentration
  int bgrids = kinetics->bgrids;

  for (int i = 0; i < bgrids; i++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuname[nu], "no2") == 0) {
        no2_nitrogen += kinetics->nus[nu][i];
      } else if (strcmp(bio->nuname[nu], "nh3") == 0) {
        nh3_nitrogen += kinetics->nus[nu][i];
      }
    }
  }

  no2_nitrogen /= bgrids;
  nh3_nitrogen /= bgrids;

  MPI_Allreduce(&no2_nitrogen,&global_no2,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&nh3_nitrogen,&global_nh3,1,MPI_DOUBLE,MPI_SUM,world);

  double diff_mass = ((global_smass - global_pre_smass) * 0.2) / (kinetics->nx * kinetics->ny * kinetics->nz);
  double diff_no2 = (global_no2 - global_pre_no2) / comm->nprocs;
  double diff_nh3 = (global_nh3 - global_pre_nh3) / comm->nprocs;

  double left = fabs(diff_mass) + fabs(diff_no2);
  double right = fabs(diff_nh3);

  if (comm->me == 0) printf("(N) Diff = %e, Biomass = %e, NO2 = %e, NH3  = %e \n",
      left-right, diff_mass, diff_no2, diff_nh3);

  global_pre_nh3 = global_nh3;
  global_pre_no2 = global_no2;
  global_pre_smass = global_smass;

  global_nh3 = 0;
  global_no2 = 0;
  global_smass = 0;
}

/* ----------------------------------------------------------------------
 biofilm benchmark problem one
 ------------------------------------------------------------------------- */

void FixVerify::benchmark_one() {
  double global_tmass, global_maxz;
  double tmass = 0;
  double maxz = 0;

  for (int i = 0; i < nlocal; i++) {
    tmass += atom->rmass[i];
    if (atom->x[i][2] > maxz) maxz = atom->x[i][2];
  }
  MPI_Allreduce(&tmass,&global_tmass,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&maxz,&global_maxz,1,MPI_DOUBLE,MPI_SUM,world);

  if (global_tmass > 1.92e-10) {
    kinetics->monod->external_gflag = 0;
  }

  if (comm->me == 0 && logfile) fprintf(logfile, "maxz = %e \n",maxz);
  if (comm->me == 0 && screen) fprintf(screen, "maxz = %e \n",maxz);

  bm1_output();
}

/* ----------------------------------------------------------------------
 biofilm benchmark problem one
 ------------------------------------------------------------------------- */

void FixVerify::benchmark_two() {
  // surface concentration
  double ssurf, gssurf, ssurf_2, ssurf_3;
  ssurf = 0;
  ssurf_2 = 0;
  ssurf_3 = 0;


  for(auto const& value: fslist) {
    int grid = kinetics->position(value);
    int up = grid + kinetics->nx * kinetics->ny;  // z direction
    int up2 = grid + (kinetics->nx * kinetics->ny)*2;  // z direction
    int up3 = grid + (kinetics->nx * kinetics->ny)*3;  // z direction
    ssurf += kinetics->nus[1][up];
    ssurf_2 += kinetics->nus[1][up2];
    ssurf_3 +=kinetics->nus[1][up3];
  }
  ssurf /= fslist.size();
  ssurf_2 /= fslist.size();
  ssurf_3 /= fslist.size();

  MPI_Allreduce(&ssurf,&gssurf,1,MPI_DOUBLE,MPI_SUM,world);

  gssurf /= comm->nprocs;
  if (comm->me == 0 && logfile) fprintf(logfile, "ssurf = %e \n",gssurf);
  if (comm->me == 0 && screen) fprintf(screen, "ssurf = %e \n",gssurf);

  if (comm->me == 0 && logfile) fprintf(logfile, "ssurf2 = %e \n",ssurf_2);
  if (comm->me == 0 && screen) fprintf(screen, "ssurf2 = %e \n",ssurf_2);

  if (comm->me == 0 && logfile) fprintf(logfile, "ssurf3 = %e \n",ssurf_3);
  if (comm->me == 0 && screen) fprintf(screen, "ssurf3 = %e \n",ssurf_3);
}


/* ----------------------------------------------------------------------
 biofilm benchmark problem three
 ------------------------------------------------------------------------- */

void FixVerify::benchmark_three() {
  double global_tmass, global_maxz;
  double tmass = 0;
  double maxz = 0;

  for (int i = 0; i < nlocal; i++) {
    tmass += atom->rmass[i];
    if (atom->x[i][2]+atom->radius[i] > maxz) maxz = atom->x[i][2]+atom->radius[i];
  }
  MPI_Allreduce(&tmass,&global_tmass,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&maxz,&global_maxz,1,MPI_DOUBLE,MPI_MAX ,world);

  if (global_maxz > 4.75e-4) {
    kinetics->monod->external_gflag = 0;
  }

  if (comm->me == 0 && logfile) fprintf(logfile, "maxz = %e \n",global_maxz);
  if (comm->me == 0 && screen) fprintf(screen, "maxz = %e \n",global_maxz);
  bm3_output();

  if (global_maxz > 5e-4) {
    kinetics->niter = -1;
    for (int i = 0; i < 50; i++) {
      kinetics->integration();
      bm3_output();
    }
    error->all(FLERR, "Finish");
  }
}

void FixVerify::free_particle_list() {
  //cutoff = 6.221e-6;
  cutoff = 1e-6;
  // uniform radius
  double d = atom->radius[0] * 2;
  double r = atom->radius[0];
  nlist.clear();
  //build neighbor list
  neighbor_list();
  // free surface particles & bottom particles
  int nfsp, gn_fsp, bp, gbp;
  double minx, miny, minz, maxx, maxy, maxz;
  double gminx, gminy, gminz, gmaxx, gmaxy, gmaxz;
  double height, ghight;

  minx = miny = minz = 10;
  maxx = maxy = 0;
  nfsp = 0;
  height = 0;
  bp = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if(nlist[i].size() == 6) continue;

    if (atom->x[i][0] < minx) minx = atom->x[i][0];
    if (atom->x[i][1] < miny) miny = atom->x[i][1];
    if (atom->x[i][2] < minz) minz = atom->x[i][2];
    if (atom->x[i][0] > maxx) maxx = atom->x[i][0];
    if (atom->x[i][1] > maxy) maxy = atom->x[i][1];

    //nfsp += 6 - nlist[i].size();
  }

  MPI_Allreduce(&minx,&gminx,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&miny,&gminy,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&minz,&gminz,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&maxx,&gmaxx,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&maxy,&gmaxy,1,MPI_DOUBLE,MPI_MAX,world);

  int gnfsurfaces;
  int nfsurfaces = 0;
  double base, top, gbase, gtop;
  base = 10;
  top = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if(nlist[i].size() == 6) continue;
    int surface = nlist[i].size();

    if (atom->x[i][0] == minx) {
      surface++;
    }
    if (atom->x[i][1] == miny) {
      surface++;
    }
    if (atom->x[i][2] == minz) {
      surface++;
      bp++;
    }
    if (atom->x[i][0] == maxx) {
      surface++;
    }
    if (atom->x[i][1] == maxy) {
      surface++;
    }

    if (surface < 6) {
      nfsp++;
      nfsurfaces += 6 - surface;
      fslist.push_back(i);
      height += atom->x[i][2];

      if (atom->x[i][2]+r > top) top = atom->x[i][2]+r;
      if (atom->x[i][2]+r < base) base = atom->x[i][2]+r;
    }
  }

  MPI_Allreduce(&nfsp,&gn_fsp,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nfsurfaces,&gnfsurfaces,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&height,&ghight,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&bp,&gbp,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&top,&gtop,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&base,&gbase,1,MPI_DOUBLE,MPI_MIN,world);

  //gn_fsp -= gnedge;

  double area = gnfsurfaces * d * d;
  double barea = gbp * d * d;

  if (comm->me == 0 && screen) fprintf(screen, "ae = %e, top = %e, base = %e, ave_h = %e \n",area/barea, gtop, gbase, ghight/(double)gn_fsp);
  if (comm->me == 0 && logfile) fprintf(logfile, "ae = %e, top = %e, base = %e, ave_h = %e \n",area/barea, gtop, gbase, ghight/(double)gn_fsp);
}


void FixVerify::neighbor_list () {
  int nall = atom->nlocal;

  for(int i = 0; i < atom->nlocal; i++){
    std::vector<int> subList;
    for(int j = 0; j < nall; j++){
      if(i != j) {
        double xd = atom->x[i][0] - atom->x[j][0];
        double yd = atom->x[i][1] - atom->x[j][1];
        double zd = atom->x[i][2] - atom->x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (atom->radius[i] + atom->radius[j] + cutoff) * (atom->radius[i] + atom->radius[j]+ cutoff);

        if (rsq <= cut) subList.push_back(j);
      }
    }
    nlist.push_back(subList);
  }
}

void FixVerify::bm1_output() {
  for (int nu = 1; nu <= nnus; nu++) {
    if (strcmp(bio->nuname[nu], "sub") == 0) {
      if (comm->me == 0 && screen) fprintf(screen, "S-sub-bulk = %e\n", kinetics->nubs[nu]);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-sub-bulk = %e\n", kinetics->nubs[nu]);
      double s = get_ave_s_sub_base();
      if (comm->me == 0 && screen) fprintf(screen, "S-sub-base = %e\n", s);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-sub-base = %e\n", s);
    }

    if (strcmp(bio->nuname[nu], "o2") == 0) {
      if (comm->me == 0 && screen) fprintf(screen, "S-o2-bulk = %e\n", kinetics->nubs[nu]);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-o2-bulk = %e\n", kinetics->nubs[nu]);
      double s = get_ave_s_o2_base();
      if (comm->me == 0 && screen) fprintf(screen, "S-o2-base = %e\n\n", s);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-o2-base = %e\n\n", s);
    }
  }
}

void FixVerify::bm3_output() {
  for (int nu = 1; nu <= nnus; nu++) {
    if (strcmp(bio->nuname[nu], "sub") == 0) {
      if (comm->me == 0 && screen) fprintf(screen, "S-sub-bulk = %e\n", kinetics->nubs[nu]);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-sub-bulk = %e\n", kinetics->nubs[nu]);
      double s = get_ave_s_sub_base();
      if (comm->me == 0 && screen) fprintf(screen, "S-sub-base = %e\n", s);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-sub-base = %e\n", s);
    }

    if (strcmp(bio->nuname[nu], "nh4") == 0) {
      if (comm->me == 0 && screen) fprintf(screen, "S-nh4-bulk = %e\n", kinetics->nubs[nu]);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-nh4-bulk = %e\n", kinetics->nubs[nu]);
      double s = get_ave_s_nh4_base();
      if (comm->me == 0 && screen) fprintf(screen, "S-nh4-base = %e\n\n", s);
      if (comm->me == 0 && logfile) fprintf(logfile, "S-nh4-base = %e\n\n", s);
    }
  }
}

double FixVerify::get_ave_s_sub_base() {
  double ave_sub_s = 0;
  double global_ave_sub_s = 0;

  int nX = kinetics->subn[0] + 2;
  int nY = kinetics->subn[1] + 2;

  for (int grid = 0; grid < kinetics->diffusion->snxx_yy_zz; grid++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuname[nu], "sub") == 0) {
        int up = grid + nX * nY;

        if (kinetics->diffusion->xgrid[grid][2] < kinetics->zlo && !kinetics->diffusion->ghost[up]) {
          ave_sub_s += kinetics->diffusion->nugrid[nu][grid];
        }
      }
    }
  }

  MPI_Allreduce(&ave_sub_s,&global_ave_sub_s,1,MPI_DOUBLE,MPI_SUM,world);

  return global_ave_sub_s/(kinetics->nx * kinetics->ny);
}

double FixVerify::get_ave_s_o2_base() {
  double ave_o2_s = 0;
  double global_ave_o2_s = 0;

  int nX = kinetics->subn[0] + 2;
  int nY = kinetics->subn[1] + 2;

  for (int grid = 0; grid < kinetics->diffusion->snxx_yy_zz; grid++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuname[nu], "o2") == 0) {
        int up = grid + nX * nY;

        if (kinetics->diffusion->xgrid[grid][2] < kinetics->zlo && !kinetics->diffusion->ghost[up]) {
          ave_o2_s += kinetics->diffusion->nugrid[nu][grid];
        }
      }
    }
  }

  MPI_Allreduce(&ave_o2_s,&global_ave_o2_s,1,MPI_DOUBLE,MPI_SUM,world);

  return global_ave_o2_s/(kinetics->nx * kinetics->ny);
}

double FixVerify::get_ave_s_nh4_base() {
  double ave_nh4_s = 0;
  double global_ave_nh4_s = 0;

  int nX = kinetics->subn[0] + 2;
  int nY = kinetics->subn[1] + 2;

  for (int grid = 0; grid < kinetics->diffusion->snxx_yy_zz; grid++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuname[nu], "nh4") == 0) {
        int up = grid + nX * nY;

        if (kinetics->diffusion->xgrid[grid][2] < kinetics->zlo && !kinetics->diffusion->ghost[up]) {
          ave_nh4_s += kinetics->diffusion->nugrid[nu][grid];
        }
      }
    }
  }

  MPI_Allreduce(&ave_nh4_s,&global_ave_nh4_s,1,MPI_DOUBLE,MPI_SUM,world);

  return global_ave_nh4_s/(kinetics->nx * kinetics->ny);
}

/* ---------------------------------------------------------------------- */

int FixVerify::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"demflag") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal fix_modify command");
    demflag = force->inumeric(FLERR, arg[1]);
    if (demflag != 0 && demflag != 1)
      error->all(FLERR, "Illegal fix divide command: demflag");
    return 2;
  }
  return 0;
}
