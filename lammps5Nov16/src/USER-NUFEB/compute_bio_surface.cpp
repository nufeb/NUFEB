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

#include "compute_bio_surface.h"

#include <math.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"
#include "domain.h"
#include "comm.h"
#include "math_const.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "stdlib.h"
#include <algorithm>

#include <vector>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

/* ---------------------------------------------------------------------- */

ComputeNufebSurface::ComputeNufebSurface(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute floc equivalent diameter command");

  scalar_flag = 1;
  extscalar = 0;

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
}

/* ---------------------------------------------------------------------- */

ComputeNufebSurface::~ComputeNufebSurface()
{
}


/* ---------------------------------------------------------------------- */
double ComputeNufebSurface::compute_scalar()
{
  if (update->ntimestep != 100000) return 0;
  invoked_scalar = update->ntimestep;
//  int ii,inum;
//  int *ilist,*numneigh;

  int localfsp = 0; //# of free surface particles
  int nfsp = 0;
  double localradi = 0;
  double nave_radi = 0;
  double ave_neighbor = 0;
  scalar = 0;
  int nneigh = 0;

  nlist.clear();

  //build neighbor list
  neighbor_list();

  for (int i = 0; i < atom->nlocal; i++) {
    nneigh += nlist[i].size();
  }

  int nnall = 0;
  MPI_Allreduce(&nneigh,&nnall,1,MPI_INT,MPI_SUM,world);

  ave_neighbor = (double)nnall/(double)atom->natoms;

  if (comm->me == 0) printf("nnal = %i, natoms = %i mean = %e \n", nnall, atom->natoms, ave_neighbor);

  for (int i = 0; i < atom->nlocal; i++) {
    localradi += atom->radius[i];
    double cutoff = atom->radius[i]*1.5;
    if (atom->x[i][0] < cutoff || atom->x[i][0] > xhi-cutoff ||
        atom->x[i][1] < cutoff || atom->x[i][1] > yhi-cutoff ||
        atom->x[i][2] < cutoff) continue;
    int jnum = nlist[i].size();
//
//    double vx, vy, vz;
//    vx = vy = vz = 0;
//
//    for (int j = 0; j < jnum; j++) {
//      double xd = atom->x[i][0] - atom->x[j][0];
//      double yd = atom->x[i][1] - atom->x[j][1];
//      double zd = atom->x[i][2] - atom->x[j][2];
//      vx += xd;
//      vy += yd;
//      vz += zd;
//    }
//    if (comm->me == 0) printf("i=%i, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e \n", i,atom->x[i][0],atom->x[i][1],atom->x[i][2],vx, vy,vz);

    if (jnum <= ave_neighbor * 29/35 ) {
      localfsp++;
      atom->type[i] = 2;
    }
  }

  MPI_Allreduce(&localradi,&nave_radi,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&localfsp,&nfsp,1,MPI_INT,MPI_SUM,world);

  if (nfsp != 0) nave_radi = nave_radi/(double)atom->natoms ;

  //if (comm->me == 0) printf("localfsp = %i, nfsp = %i \n", localfsp, nfsp);

  scalar = (3.1415 * nave_radi * nave_radi * nfsp) / (yhi * xhi);
  //printf("nfsp = %i sa = %e ae = %e \n", nfsp, scalar, scalar * (yhi * xhi));

  return scalar;
}

void ComputeNufebSurface::neighbor_list () {
  int nall = atom->nlocal + atom->nghost;

  for(int i = 0; i < atom->nlocal; i++){
    std::vector<int> subList;
    for(int j = 0; j < nall; j++){
      if(i != j) {
        double xd = atom->x[i][0] - atom->x[j][0];
        double yd = atom->x[i][1] - atom->x[j][1];
        double zd = atom->x[i][2] - atom->x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (atom->radius[i] + atom->radius[j] + 1.0e-6) * (atom->radius[i] + atom->radius[j]+ 1.0e-6);

        if (rsq <= cut) subList.push_back(j);
      }
    }
    nlist.push_back(subList);
  }
}
