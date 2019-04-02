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

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include "mpi.h"
#include "stdio.h"
#include <vector>
#include <algorithm>
/* ifdef allow this file to be included in a C program */

#ifdef ENABLE_CFDDEM
#ifdef __cplusplus
extern "C" {
#endif

  void lammps_open(int, char **, MPI_Comm, void **);  /* start-up LAMMPS */
  void lammps_close(void *);                          /* shut-down LAMMPS */
  void lammps_file(void *, char *);                 /* run input script */
  char *lammps_command(void *, char *);             /* execute a command */

  void lammps_sync(void* ptr);
  int lammps_get_global_n(void *);              /* return # of atoms */

  /* get atom x&v&d & rho & type */
  void lammps_get_initial_np(void* ptr, int* np_);

  void lammps_get_initial_info(void* ptr, double* coords, double* velos,
                               double* diam, double* rho_, int* tag_,
                               int* lmpCpuId_, int* type_);

  /* get atom number in each lmp cpu */
  int lammps_get_local_n(void* ptr);

  /* get local domain in each lmp cpu */
  void lammps_get_local_domain(void* ptr, double* domain_);

  /* get atom x&v&foamCpuId&lmpCpuId */
  void lammps_get_local_info(void* ptr, double* coords_, double* velos_,
                             int* foamCpuId_, int* lmpCpuId_, int* tag_,
  						   double* rho_, double* diam_, int* type_);

  /* set atom x&v&foamCpuId for all procs */
  void lammps_put_local_info(void* ptr, int nLocalIn, double* fdrag, 
                             double* DuDt, int* foamCpuIdIn,
							 int* tagIn, double* gridPU);

  void lammps_step(void* ptr, int n);
  void lammps_step_stop(void* ptr, int n, int nsteps);
  void lammps_set_timestep(void* ptr, double dt_i);
  double lammps_get_timestep(void* ptr);
  int lammps_get_demflag(void *ptr);
  void lammps_set_demflag(void *ptr, int x);

  int lammps_get_bio_steps(void* ptr);
  int lammps_get_dem_steps(void* ptr);
  int lammps_get_nloops(void* ptr);
  double lammps_get_bio_dt(void* ptr);
  double lammps_get_dem_dt(void* ptr);
  double lammps_get_scale_time(void* ptr);

  void lammps_create_particle(void* ptr, int npAdd, double* position, double* tag, 
                              double diameter, double rho, int type, double* vel);
  void lammps_delete_particle(void* ptr, int* deleteList, int nDelete);

  /* used in the sorting part when assigning data from OpenFOAM*/
  struct tagpair {
    int tag;
    int index;
  };

  struct by_number {
    bool operator() (tagpair const &left, tagpair const &right) {
      return left.tag < right.tag;
    }
  };
#ifdef __cplusplus
}
#endif
#endif

