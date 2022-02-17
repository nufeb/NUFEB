/***************************************************************************
                                   dpd.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the dpd pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Jan 15, 2014
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "dpd_cl.h"
#elif defined(USE_CUDART)
const char *dpd=0;
#else
#include "dpd_cubin.h"
#endif

#include "lal_dpd.h"
#include <cassert>
namespace LAMMPS_AL {
#define DPDT DPD<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
DPDT::DPD() : BaseDPD<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
DPDT::~DPD() {
  clear();
}

template <class numtyp, class acctyp>
int DPDT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int DPDT::init(const int ntypes,
               double **host_cutsq, double **host_a0,
               double **host_gamma, double **host_sigma,
               double **host_cut, double *host_special_lj,
               const bool tstat_only,
               const int nlocal, const int nall,
               const int max_nbors, const int maxspecial,
               const double cell_size,
               const double gpu_split, FILE *_screen) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,_screen,dpd,"k_dpd");
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  _lj_types=lj_types;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  coeff.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,coeff,host_write,host_a0,host_gamma,
                         host_sigma,host_cut);

  UCL_H_Vec<numtyp> host_rsq(lj_types*lj_types,*(this->ucl_device),
                             UCL_WRITE_ONLY);
  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_rsq,host_cutsq);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _tstat_only = 0;
  if (tstat_only) _tstat_only=1;

  _allocated=true;
  this->_max_bytes=coeff.row_bytes()+cutsq.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void DPDT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  coeff.clear();
  cutsq.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double DPDT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(DPD<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void DPDT::loop(const bool _eflag, const bool _vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;

  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_fast.set_size(GX,BX);
    this->k_pair_fast.run(&this->atom->x, &coeff, &sp_lj,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->v, &cutsq,
                          &this->_dtinvsqrt, &this->_seed, &this->_timestep,
                          &this->_tstat_only, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &coeff, &_lj_types, &sp_lj,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->atom->v, &cutsq, &this->_dtinvsqrt,
                     &this->_seed, &this->_timestep, &this->_tstat_only,
                     &this->_threads_per_atom);
  }
  this->time_pair.stop();
}

template <class numtyp, class acctyp>
void DPDT::update_coeff(int ntypes, double **host_a0, double **host_gamma,
                        double **host_sigma, double **host_cut)
{
  UCL_H_Vec<numtyp> host_write(_lj_types*_lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);
  this->atom->type_pack4(ntypes,_lj_types,coeff,host_write,host_a0,host_gamma,
                         host_sigma,host_cut);
}

template class DPD<PRECISION,ACC_PRECISION>;
}
