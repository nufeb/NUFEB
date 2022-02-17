/***************************************************************************
                                  atom.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for particle data management

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_atom.h"

namespace LAMMPS_AL {
#define AtomT Atom<numtyp,acctyp>

template <class numtyp, class acctyp>
AtomT::Atom() : _compiled(false),_allocated(false),
                              _max_gpu_bytes(0) {
  #ifdef USE_CUDPP
  sort_config.op = CUDPP_ADD;
  sort_config.datatype = CUDPP_UINT;
  sort_config.algorithm = CUDPP_SORT_RADIX;
  sort_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
  #endif
}

template <class numtyp, class acctyp>
int AtomT::bytes_per_atom() const {
  int id_space=0;
  if (_gpu_nbor==1)
    id_space=2;
  else if (_gpu_nbor==2)
    id_space=4;
  int bytes=4*sizeof(numtyp)+id_space*sizeof(int);
  if (_rot)
    bytes+=4*sizeof(numtyp);
  if (_charge)
    bytes+=sizeof(numtyp);
  if (_vel)
    bytes+=4*sizeof(numtyp);
  return bytes;
}

template <class numtyp, class acctyp>
bool AtomT::alloc(const int nall) {
  _max_atoms=static_cast<int>(static_cast<double>(nall)*1.10);

  bool success=true;

  // Ignore host/device transfers?
  _host_view=false;
  if (dev->shared_memory() && sizeof(numtyp)==sizeof(double)) {
    _host_view=true;
    #ifdef GPU_CAST
    assert(0==1);
    #endif
  }

  // Allocate storage for CUDPP sort
  #ifdef USE_CUDPP
  if (_gpu_nbor==1) {
    CUDPPResult result = cudppPlan(&sort_plan, sort_config, _max_atoms, 1, 0);
    if (CUDPP_SUCCESS != result)
      return false;
  }
  #endif

  // ---------------------------  Device allocations
  int gpu_bytes=0;
  success=success && (x.alloc(_max_atoms*4,*dev,UCL_WRITE_ONLY,
                              UCL_READ_ONLY)==UCL_SUCCESS);
  #ifdef GPU_CAST
  success=success && (x_cast.alloc(_max_atoms*3,*dev,UCL_WRITE_ONLY,
                                   UCL_READ_ONLY)==UCL_SUCCESS);
  success=success && (type_cast.alloc(_max_atoms,*dev,UCL_WRITE_ONLY,
                                      UCL_READ_ONLY)==UCL_SUCCESS);
  gpu_bytes+=x_cast.device.row_bytes()+type_cast.device.row_bytes();
  #endif

  if (_charge && _host_view==false) {
    success=success && (q.alloc(_max_atoms,*dev,UCL_WRITE_ONLY,
                                UCL_READ_ONLY)==UCL_SUCCESS);
    gpu_bytes+=q.device.row_bytes();
  }
  if (_rot && _host_view==false) {
    success=success && (quat.alloc(_max_atoms*4,*dev,UCL_WRITE_ONLY,
                                   UCL_READ_ONLY)==UCL_SUCCESS);
    gpu_bytes+=quat.device.row_bytes();
  }
  if (_vel && _host_view==false) {
    success=success && (v.alloc(_max_atoms*4,*dev,UCL_WRITE_ONLY,
                                   UCL_READ_ONLY)==UCL_SUCCESS);
    gpu_bytes+=v.device.row_bytes();
  }

  if (_gpu_nbor>0) {
    if (_bonds) {
      success=success && (dev_tag.alloc(_max_atoms,*dev,
                                        UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=dev_tag.row_bytes();
    }
    if (_gpu_nbor==1) {
      success=success && (dev_cell_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
      gpu_bytes+=dev_cell_id.row_bytes();
    } else {
      success=success && (host_particle_id.alloc(_max_atoms,*dev,
                                                 UCL_WRITE_ONLY)==UCL_SUCCESS);
      success=success &&
             (host_cell_id.alloc(_max_atoms,*dev,UCL_NOT_PINNED)==UCL_SUCCESS);
    }
    if (_gpu_nbor==2 && _host_view)
      dev_particle_id.view(host_particle_id);
    else
      success=success && (dev_particle_id.alloc(_max_atoms,*dev,
                                                UCL_READ_ONLY)==UCL_SUCCESS);
    gpu_bytes+=dev_particle_id.row_bytes();
  }

  gpu_bytes+=x.device.row_bytes();
  if (gpu_bytes>_max_gpu_bytes)
    _max_gpu_bytes=gpu_bytes;

  _allocated=true;
  return success;
}

template <class numtyp, class acctyp>
bool AtomT::add_fields(const bool charge, const bool rot,
                       const int gpu_nbor, const bool bonds, const bool vel) {
  bool success=true;
  // Ignore host/device transfers?
  int gpu_bytes=0;

  if (charge && _charge==false) {
    _charge=true;
    _other=true;
    if (_host_view==false) {
      success=success && (q.alloc(_max_atoms,*dev,UCL_WRITE_ONLY,
                                  UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=q.device.row_bytes();
    }
  }

  if (rot && _rot==false) {
    _rot=true;
    _other=true;
    if (_host_view==false) {
      success=success && (quat.alloc(_max_atoms*4,*dev,UCL_WRITE_ONLY,
                                     UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=quat.device.row_bytes();
    }
  }

  if (vel && _vel==false) {
    _vel=true;
    _other=true;
    if (_host_view==false) {
      success=success && (v.alloc(_max_atoms*4,*dev,UCL_WRITE_ONLY,
                                     UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=v.device.row_bytes();
    }
  }

  if (bonds && _bonds==false) {
    _bonds=true;
    if (_bonds && _gpu_nbor>0) {
      success=success && (dev_tag.alloc(_max_atoms,*dev,
                                        UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=dev_tag.row_bytes();
    }
  }

  if (gpu_nbor>0 && _gpu_nbor==0) {
    _gpu_nbor=gpu_nbor;
    #ifdef USE_CUDPP
    if (_gpu_nbor==1) {
      CUDPPResult result = cudppPlan(&sort_plan, sort_config, _max_atoms, 1, 0);
      if (CUDPP_SUCCESS != result)
        return false;
    }
    #endif
    success=success && (dev_particle_id.alloc(_max_atoms,*dev,
                                              UCL_READ_ONLY)==UCL_SUCCESS);
    gpu_bytes+=dev_particle_id.row_bytes();
    if (_bonds) {
      success=success && (dev_tag.alloc(_max_atoms,*dev,
                                        UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=dev_tag.row_bytes();
    }
    if (_gpu_nbor==1) {
      success=success && (dev_cell_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
      gpu_bytes+=dev_cell_id.row_bytes();
    } else {
      success=success && (host_particle_id.alloc(_max_atoms,*dev,
                                                 UCL_WRITE_ONLY)==UCL_SUCCESS);
      success=success &&
             (host_cell_id.alloc(_max_atoms,*dev,UCL_NOT_PINNED)==UCL_SUCCESS);
    }
  }

  return success;
}

template <class numtyp, class acctyp>
bool AtomT::init(const int nall, const bool charge, const bool rot,
                 UCL_Device &devi, const int gpu_nbor, const bool bonds, const bool vel) {
  clear();

  bool success=true;
  _x_avail=false;
  _q_avail=false;
  _quat_avail=false;
  _v_avail=false;
  _resized=false;
  _gpu_nbor=gpu_nbor;
  _bonds=bonds;
  _charge=charge;
  _rot=rot;
  _vel=vel;
  _other=_charge || _rot || _vel;
  dev=&devi;
  _time_transfer=0;

  // Initialize atom and nbor data
  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;

  // Initialize timers for the selected device
  time_pos.init(*dev);
  time_q.init(*dev);
  time_quat.init(*dev);
  time_vel.init(*dev);
  time_pos.zero();
  time_q.zero();
  time_quat.zero();
  time_vel.zero();
  _time_cast=0.0;

  #ifdef GPU_CAST
  compile_kernels(*dev);
  #endif

  return success && alloc(ef_nall);
}

template <class numtyp, class acctyp>
void AtomT::clear_resize() {
  if (!_allocated)
    return;
  _allocated=false;

  x.clear();
  if (_charge)
    q.clear();
  if (_rot)
    quat.clear();
  if (_vel)
    v.clear();

  dev_cell_id.clear();
  dev_particle_id.clear();
  dev_tag.clear();
  #ifdef GPU_CAST
  x_cast.clear();
  type_cast.clear();
  #endif

  #ifdef USE_CUDPP
  if (_gpu_nbor==1) cudppDestroyPlan(sort_plan);
  #endif

  if (_gpu_nbor==2) {
    host_particle_id.clear();
    host_cell_id.clear();
  }
}

template <class numtyp, class acctyp>
void AtomT::clear() {
  _max_gpu_bytes=0;
  if (!_allocated)
    return;

  time_pos.clear();
  time_q.clear();
  time_quat.clear();
  time_vel.clear();
  clear_resize();

  #ifdef GPU_CAST
  if (_compiled) {
    k_cast_x.clear();
    delete atom_program;
    _compiled=false;
  }
  #endif
}

template <class numtyp, class acctyp>
double AtomT::host_memory_usage() const {
  int atom_bytes=4;
  if (_charge)
    atom_bytes+=1;
  if (_rot)
    atom_bytes+=4;
  if (_vel)
    atom_bytes+=4;
  return _max_atoms*atom_bytes*sizeof(numtyp)+sizeof(Atom<numtyp,acctyp>);
}

// Sort arrays for neighbor list calculation
template <class numtyp, class acctyp>
void AtomT::sort_neighbor(const int num_atoms) {
  #ifdef USE_CUDPP
  CUDPPResult result = cudppSort(sort_plan, (unsigned *)dev_cell_id.begin(),
                                 (int *)dev_particle_id.begin(),
                                 8*sizeof(unsigned), num_atoms);
  if (CUDPP_SUCCESS != result) {
    printf("Error in cudppSort\n");
    UCL_GERYON_EXIT;
  }
  #endif
}

#ifdef GPU_CAST
#if defined(USE_OPENCL)
#include "atom_cl.h"
#elif defined(USE_CUDART)
const char *atom=0;
#else
#include "atom_cubin.h"
#endif

template <class numtyp, class acctyp>
void AtomT::compile_kernels(UCL_Device &dev) {
  std::string flags = "-D"+std::string(OCL_VENDOR);
  atom_program=new UCL_Program(dev);
  atom_program->load_string(atom,flags);
  k_cast_x.set_function(*atom_program,"kernel_cast_x");
  _compiled=true;
}

#endif

template class Atom<PRECISION,ACC_PRECISION>;
}
