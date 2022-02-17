/***************************************************************************
                                 mie_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to mie acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>

#include "lal_mie.h"

using namespace std;
using namespace LAMMPS_AL;

static Mie<PRECISION,ACC_PRECISION> MLMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int mie_gpu_init(const int ntypes, double **cutsq, double **host_mie1,
                 double **host_mie2, double **host_mie3, double **host_mie4,
                 double **host_gamA, double **host_gamR,
                 double **offset, double *special_lj,
                 const int inum, const int nall, const int max_nbors,
                 const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen) {
  MLMF.clear();
  gpu_mode=MLMF.device->gpu_mode();
  double gpu_split=MLMF.device->particle_split();
  int first_gpu=MLMF.device->first_device();
  int last_gpu=MLMF.device->last_device();
  int world_me=MLMF.device->world_me();
  int gpu_rank=MLMF.device->gpu_rank();
  int procs_per_gpu=MLMF.device->procs_per_gpu();

  MLMF.device->init_message(screen,"mie",first_gpu,last_gpu);

  bool message=false;
  if (MLMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=MLMF.init(ntypes, cutsq, host_mie1, host_mie2,
                      host_mie3, host_mie4, host_gamA, host_gamR,
                      offset, special_lj, inum, nall, 300,
                      maxspecial, cell_size, gpu_split, screen);

  MLMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing Device %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing Devices %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      init_ok=MLMF.init(ntypes, cutsq, host_mie1, host_mie2,
                        host_mie3, host_mie4, host_gamA, host_gamR,
                        offset, special_lj, inum, nall, 300, maxspecial,
                        cell_size, gpu_split, screen);

    MLMF.device->gpu_barrier();
    if (message)
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    MLMF.estimate_gpu_overhead();
  return init_ok;
}

void mie_gpu_clear() {
  MLMF.clear();
}

int ** mie_gpu_compute_n(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success) {
  return MLMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                      subhi, tag, nspecial, special, eflag, vflag, eatom,
                      vatom, host_start, ilist, jnum, cpu_time, success);
}

void mie_gpu_compute(const int ago, const int inum_full, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success) {
  MLMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,
               firstneigh,eflag,vflag,eatom,vatom,host_start,cpu_time,success);
}

double mie_gpu_bytes() {
  return MLMF.host_memory_usage();
}


