/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_bio_utilities.h"

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "mpi.h"
#include "comm.h"
#include "memory.h"
#include "input.h"
#include "variable.h"
#include <iostream>
#include "math_const.h"
#include <vector>
#include <algorithm>


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixUtilities::FixUtilities(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix utilities command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix utilities command: calling steps should be positive integer");
}

FixUtilities::~FixUtilities()
{
}

/* ---------------------------------------------------------------------- */

int FixUtilities::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUtilities::init()
{
}


void FixUtilities::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  nall = atom->natoms;

  int* tags;
  double* r;
  double* flatten_x;
  double** x;

  if(comm->me == 0){
    tags = memory->create(tags, nall, "utility/all_tags");
    r = memory->create(r, nall, "utility/all_tags");
    flatten_x = memory->create(flatten_x, nall*3, "utility/all_x");
  }

  // pass all required info to proc 0
  comunicator(tags, flatten_x, r);

  if (comm->me == 0) {
    x = memory->create(x, nall, 3, "utility/x");

    for (int i = 0; i < nall; i++){
      x[i][0] = flatten_x[i*3];
      x[i][1] = flatten_x[i*3+1];
      x[i][2] = flatten_x[i*3+2];
     // printf("i=%i, x=%e, y=%e, z=%e, r=%e, tag=%i \n", i, flatten_x[i*3], flatten_x[i*3+1], flatten_x[i*3+2], r[i], tags[i]);
    }

    pFile = fopen ("floc", "a");

    list.clear();
    fourThirdsPI = 4.0*MY_PI/3.0;

    visit = new int[nall]();
    vector <vector<int>> flocs;
    vector<int> parent_floc;
    int max_size = 0;

    //build neighbor list
    neighbor_list(x, r);

    for (int i = 0; i < nall; i++) {
      if (visit[i] == 0) {
	vector<int> floc;
	floc.push_back(i);
	get_floc (i, floc);

	//get parent floc,assuming floc with largest size is parent
	if (floc.size() > max_size) {
	  max_size = floc.size();
	  parent_floc = floc;
	}
	flocs.push_back(floc);
      }
    }

    fprintf(pFile, "ntimestep = %i \n", update->ntimestep);
    fprintf(pFile, "Parent floc size = %i \n", parent_floc.size());
    //remove parent floc
    flocs.erase(remove( flocs.begin(), flocs.end(), parent_floc ), flocs.end() );

    for (const vector<int> &floc : flocs) {
      //new detached floc
      int new_size = 0;
      int old_size = 0;

      double old_vol = 0;
      double new_vol = 0;
      double floc_x, floc_y ,floc_z;

      for (const int &i : floc) {
	double radius = r[i];
	double volume = fourThirdsPI * radius * radius * radius;
	//i not in id
	if (find(floc_tags.begin(), floc_tags.end(), tags[i]) == floc_tags.end()) {
	  new_size++;
	  floc_tags.push_back(tags[i]);
	  new_vol = volume + new_vol;
	} else {
	  old_size++;
	  old_vol = volume + old_vol;
	}
	floc_x += x[i][0];
	floc_y += x[i][1];
	floc_z += x[i][2];
      }

      int size = floc.size();
      floc_x = floc_x/size;
      floc_y = floc_y/size;
      floc_z = floc_z/size;

      if (old_size == size) continue;
      else if (new_size == size) {
	fprintf(pFile, "New floc at ntimestep = %i \n", update->ntimestep);
	fprintf(pFile, "  Average position x = %e, y = %e, z = %e \n", floc_x, floc_y, floc_z);
	fprintf(pFile, "  Size = %i \n", size);
	fprintf(pFile, "  Volume = %e \n", new_vol);
	fprintf(pFile, "\n");
      } else {
	fprintf(pFile, "Mixed floc at ntimestep = %i \n", update->ntimestep);
	fprintf(pFile, "  Average position x = %e, y = %e, z = %e \n", floc_x, floc_y, floc_z);
	fprintf(pFile, "  Size = %i \n", size);
	fprintf(pFile, "  Old size = %i \n", old_size);
	fprintf(pFile, "  New size = %i \n", new_size);
	fprintf(pFile, "  Volume = %e \n", new_vol + old_vol);
	fprintf(pFile, "  Old volume = %e \n", old_vol);
	fprintf(pFile, "  New volume = %e \n", new_vol);
	fprintf(pFile, "\n");
      }
    }

    fclose(pFile);

    memory->destroy(x);
    memory->destroy(tags);
    memory->destroy(r);
    memory->destroy(flatten_x);
    delete[] visit;
  }
}

/* ----------------------------------------------------------------------
 Gathering all required information in Proc 0
 ------------------------------------------------------------------------- */
void FixUtilities::comunicator (int* tags, double* flatten_x, double* r) {
  int nlocal = atom->nlocal;

  double* flatten_local_x;

  flatten_local_x = memory->create(flatten_local_x, nlocal*3, "utility/all_x");

  for (int i = 0; i < nlocal; i++){
    flatten_local_x[i*3] = atom->x[i][0];
    flatten_local_x[i*3+1] = atom->x[i][1];
    flatten_local_x[i*3+2] = atom->x[i][2];
  }

  int *recvcounts = new int[comm->nprocs];
  int *recvcounts2 = new int[comm->nprocs];
  int *displs = new int[comm->nprocs];
  int *displs2 = new int[comm->nprocs];
  // setup for gatherv
  int nlocal2 = nlocal*3;
  MPI_Allgather(&nlocal,1,MPI_INT,recvcounts,1,MPI_INT,world);
  MPI_Allgather(&nlocal2,1,MPI_INT,recvcounts2,1,MPI_INT,world);

  displs[0] = 0;
  displs2[0] = 0;

  for (int iproc = 1; iproc < comm->nprocs; iproc++) {
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
    displs2[iproc] = displs2[iproc-1] + recvcounts2[iproc-1];
  }

  MPI_Gatherv(atom->tag, nlocal, MPI_INT, tags, recvcounts, displs, MPI_INT, 0, world);
  MPI_Gatherv(atom->radius, nlocal, MPI_DOUBLE, r, recvcounts, displs, MPI_DOUBLE, 0, world);
  MPI_Gatherv(flatten_local_x, nlocal*3, MPI_DOUBLE, flatten_x, recvcounts2, displs2, MPI_DOUBLE, 0, world);

  memory->destroy(flatten_local_x);
}

/* ----------------------------------------------------------------------*/
void FixUtilities::get_floc (int bac, vector<int>& floc) {
  visit[bac] = 1;

  for (int const& j: list.at(bac)) {

    if (visit[j] == 0) {
      floc.push_back(j);
      get_floc (j, floc);
    }
  }
}

/* ----------------------------------------------------------------------*/
void FixUtilities::neighbor_list (double** x, double* r) {

  for(int i = 0; i < nall; i++){
    vector<int> subList;
    for(int j = 0; j < nall; j++){
      if(i != j) {
        double xd = x[i][0] - x[j][0];
        double yd = x[i][1] - x[j][1];
        double zd = x[i][2] - x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (r[i] + r[j] + 1.0e-6) * (r[i] + r[j]+ 1.0e-6);

        if (rsq <= cut) subList.push_back(j);
      }
    }
    list.push_back(subList);
  }
}
