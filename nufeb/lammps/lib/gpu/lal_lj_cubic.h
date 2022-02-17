/***************************************************************************
                                 lj_cubic.h
                             -------------------
                              Trung Dac Nguyen

  Class for acceleration of the lj/cubic pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#ifndef LAL_LJ_H
#define LAL_LJ_H

#include "lal_base_atomic.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class LJCubic : public BaseAtomic<numtyp, acctyp> {
 public:
  LJCubic();
  ~LJCubic();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    *
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init(const int ntypes, double **host_cutsq, double **host_cut_inner_sq,
           double **host_cut_inner, double **host_sigma, double **host_epsilon,
           double **host_lj1, double **host_lj2, double **host_lj3,
           double **host_lj4, double *host_special_lj, const int nlocal,
           const int nall, const int max_nbors, const int maxspecial,
           const double cell_size, const double gpu_split, FILE *screen);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// lj1.x = lj1, lj1.y = lj2, lj1.z = cutsq
  UCL_D_Vec<numtyp4> lj1;
  /// lj2.x = cut_inner_sq, lj2.y = cut_inner, lj2.z = sigma, lj2.w = epsilon
  UCL_D_Vec<numtyp4> lj2;
  /// lj3.x = lj3, lj3.y = lj4
  UCL_D_Vec<numtyp2> lj3;
  /// Special LJ values
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

 private:
  bool _allocated;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
