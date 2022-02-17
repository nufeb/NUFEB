/***************************************************************************
                                 tersoff.cpp
                             -------------------
                               Trung Dac Nguyen

  Class for acceleration of the tersoff pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Thu April 17, 2014
    email                : ndactrung@gmail.com
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "tersoff_cl.h"
#elif defined(USE_CUDART)
const char *tersoff=0;
#else
#include "tersoff_cubin.h"
#endif

#include "lal_tersoff.h"
#include <cassert>
namespace LAMMPS_AL {
#define TersoffT Tersoff<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
TersoffT::Tersoff() : BaseThree<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
TersoffT::~Tersoff() {
  clear();
}

template <class numtyp, class acctyp>
int TersoffT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int TersoffT::init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
                   const double cell_size, const double gpu_split, FILE *_screen,
                   int* host_map, const int nelements, int*** host_elem2param, const int nparams,
                   const double* lam1, const double* lam2, const double* lam3,const double* powermint,
                   const double* biga, const double* bigb, const double* bigr, const double* bigd,
                   const double* c1, const double* c2, const double* c3, const double* c4,
                   const double* c, const double* d, const double* h, const double* gamma,
                   const double* beta, const double* powern, const double* host_cutsq)
{
  int success;
  success=this->init_three(nlocal,nall,max_nbors,0,cell_size,gpu_split,
                           _screen,tersoff,"k_tersoff_repulsive",
                           "k_tersoff_three_center", "k_tersoff_three_end",
                           "k_tersoff_short_nbor");
  if (success!=0)
    return success;

  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;
  _zetaij.alloc(ef_nall*max_nbors,*(this->ucl_device),UCL_READ_WRITE);

  k_zeta.set_function(*(this->pair_program),"k_tersoff_zeta");

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  _lj_types=lj_types;

  _nparams = nparams;
  _nelements = nelements;

  UCL_H_Vec<numtyp4> dview(nparams,*(this->ucl_device),
                           UCL_WRITE_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=(numtyp)0;
    dview[i].y=(numtyp)0;
    dview[i].z=(numtyp)0;
    dview[i].w=(numtyp)0;
  }

  // pack coefficients into arrays
  ts1.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(lam1[i]);
    dview[i].y=static_cast<numtyp>(lam2[i]);
    dview[i].z=static_cast<numtyp>(lam3[i]);
    dview[i].w=static_cast<numtyp>(powermint[i]);
  }

  ucl_copy(ts1,dview,false);
  ts1_tex.get_texture(*(this->pair_program),"ts1_tex");
  ts1_tex.bind_float(ts1,4);

  ts2.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(biga[i]);
    dview[i].y=static_cast<numtyp>(bigb[i]);
    dview[i].z=static_cast<numtyp>(bigr[i]);
    dview[i].w=static_cast<numtyp>(bigd[i]);
  }

  ucl_copy(ts2,dview,false);
  ts2_tex.get_texture(*(this->pair_program),"ts2_tex");
  ts2_tex.bind_float(ts2,4);

  ts3.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(c1[i]);
    dview[i].y=static_cast<numtyp>(c2[i]);
    dview[i].z=static_cast<numtyp>(c3[i]);
    dview[i].w=static_cast<numtyp>(c4[i]);
  }

  ucl_copy(ts3,dview,false);
  ts3_tex.get_texture(*(this->pair_program),"ts3_tex");
  ts3_tex.bind_float(ts3,4);

  ts4.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(c[i]);
    dview[i].y=static_cast<numtyp>(d[i]);
    dview[i].z=static_cast<numtyp>(h[i]);
    dview[i].w=static_cast<numtyp>(gamma[i]);
  }

  ucl_copy(ts4,dview,false);
  ts4_tex.get_texture(*(this->pair_program),"ts4_tex");
  ts4_tex.bind_float(ts4,4);

  ts5.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(beta[i]);
    dview[i].y=static_cast<numtyp>(powern[i]);
    dview[i].z=(numtyp)0;
    dview[i].w=(numtyp)0;
  }

  ucl_copy(ts5,dview,false);
  ts5_tex.get_texture(*(this->pair_program),"ts5_tex");
  ts5_tex.bind_float(ts5,4);

  UCL_H_Vec<numtyp> cutsq_view(nparams,*(this->ucl_device),
                               UCL_WRITE_ONLY);
  double cutsqmax = 0.0;
  for (int i=0; i<nparams; i++) {
    cutsq_view[i]=static_cast<numtyp>(host_cutsq[i]);
    if (cutsqmax < host_cutsq[i]) cutsqmax = host_cutsq[i];
  }
  cutsq.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(cutsq,cutsq_view,false);

  _cutshortsq = static_cast<numtyp>(cutsqmax);

  UCL_H_Vec<int> dview_elem2param(nelements*nelements*nelements,
                           *(this->ucl_device), UCL_WRITE_ONLY);

  elem2param.alloc(nelements*nelements*nelements,*(this->ucl_device),
                   UCL_READ_ONLY);

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++) {
         int idx = i*nelements*nelements+j*nelements+k;
         dview_elem2param[idx] = host_elem2param[i][j][k];
      }

  ucl_copy(elem2param,dview_elem2param,false);

  UCL_H_Vec<int> dview_map(lj_types, *(this->ucl_device), UCL_WRITE_ONLY);
  for (int i = 0; i < ntypes; i++)
    dview_map[i] = host_map[i];

  map.alloc(lj_types,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(map,dview_map,false);

  _allocated=true;
  this->_max_bytes=ts1.row_bytes()+ts2.row_bytes()+ts3.row_bytes()+
    ts4.row_bytes()+ts5.row_bytes()+cutsq.row_bytes()+
    map.row_bytes()+elem2param.row_bytes()+_zetaij.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void TersoffT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  ts1.clear();
  ts2.clear();
  ts3.clear();
  ts4.clear();
  ts5.clear();
  cutsq.clear();
  map.clear();
  elem2param.clear();
  _zetaij.clear();

  k_zeta.clear();

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double TersoffT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Tersoff<numtyp,acctyp>);
}

#define KTHREADS this->_threads_per_atom
#define JTHREADS this->_threads_per_atom
// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void TersoffT::loop(const bool _eflag, const bool _vflag, const int evatom) {
  // Compute the block size and grid size to keep all cores busy
  int BX=this->block_pair();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;

  // build the short neighbor list
  int ainum=this->_ainum;
  int nbor_pitch=this->nbor->nbor_pitch();
  int GX=static_cast<int>(ceil(static_cast<double>(ainum)/
                               (BX/this->_threads_per_atom)));

  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &cutsq, &map,
                 &elem2param, &_nelements, &_nparams,
                 &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                 &this->dev_short_nbor, &ainum,
                 &nbor_pitch, &this->_threads_per_atom);

  // re-allocate zetaij if necessary
  int nall = this->_nall;
  if (nall*this->_max_nbors > _zetaij.cols()) {
    int _nmax=static_cast<int>(static_cast<double>(nall)*1.10);
    _zetaij.resize(this->_max_nbors*_nmax);
  }

  nbor_pitch=this->nbor->nbor_pitch();
  GX=static_cast<int>(ceil(static_cast<double>(this->_ainum)/
                               (BX/(JTHREADS*KTHREADS))));

  this->k_zeta.set_size(GX,BX);
  this->k_zeta.run(&this->atom->x, &ts1, &ts2, &ts3, &ts4, &ts5, &cutsq,
                   &map, &elem2param, &_nelements, &_nparams, &_zetaij,
                   &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                   &this->dev_short_nbor,
                   &eflag, &this->_ainum, &nbor_pitch, &this->_threads_per_atom);

  ainum=this->ans->inum();
  nbor_pitch=this->nbor->nbor_pitch();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  this->time_pair.start();
  this->k_pair.set_size(GX,BX);
  this->k_pair.run(&this->atom->x, &ts1, &ts2, &cutsq,
                   &map, &elem2param, &_nelements, &_nparams,
                   &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                   &this->dev_short_nbor,
                   &this->ans->force, &this->ans->engv,
                   &eflag, &vflag, &ainum, &nbor_pitch,
                   &this->_threads_per_atom);

  BX=this->block_size();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                           (BX/(KTHREADS*JTHREADS))));
  this->k_three_center.set_size(GX,BX);
  this->k_three_center.run(&this->atom->x, &ts1, &ts2, &ts4, &cutsq,
                           &map, &elem2param, &_nelements, &_nparams, &_zetaij,
                           &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                           &this->dev_short_nbor,
                           &this->ans->force, &this->ans->engv, &eflag, &vflag, &ainum,
                           &nbor_pitch, &this->_threads_per_atom, &evatom);

  Answer<numtyp,acctyp> *end_ans;
  #ifdef THREE_CONCURRENT
  end_ans=this->ans2;
  #else
  end_ans=this->ans;
  #endif
  if (evatom!=0) {
    this->k_three_end_vatom.set_size(GX,BX);
    this->k_three_end_vatom.run(&this->atom->x, &ts1, &ts2, &ts4, &cutsq,
                          &map, &elem2param, &_nelements, &_nparams, &_zetaij,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->nbor->dev_ilist, &this->dev_short_nbor,
                          &end_ans->force, &end_ans->engv, &eflag, &vflag, &ainum,
                          &nbor_pitch, &this->_threads_per_atom, &this->_gpu_nbor);

  } else {
    this->k_three_end.set_size(GX,BX);
    this->k_three_end.run(&this->atom->x, &ts1, &ts2, &ts4, &cutsq,
                          &map, &elem2param, &_nelements, &_nparams, &_zetaij,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->nbor->dev_ilist, &this->dev_short_nbor,
                          &end_ans->force, &end_ans->engv, &eflag, &vflag, &ainum,
                          &nbor_pitch, &this->_threads_per_atom, &this->_gpu_nbor);
  }

  this->time_pair.stop();
}

template class Tersoff<PRECISION,ACC_PRECISION>;
}
