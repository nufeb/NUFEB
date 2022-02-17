/***************************************************************************
                                  coul_dsf.h
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the coul/dsf pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 8/15/2012
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#ifndef LAL_LJ_DSF_H
#define LAL_LJ_DSF_H

#include "lal_base_charge.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class CoulDSF : public BaseCharge<numtyp, acctyp> {
 public:
  CoulDSF();
  ~CoulDSF();

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
  int init(const int ntypes, const int nlocal, const int nall,
           const int max_nbors, const int maxspecial,
           const double cell_size, const double gpu_split, FILE *screen,
           const double host_cut_coulsq, double *host_special_coul,
           const double qqrd2e, const double e_shift, const double f_shift,
           const double alpha);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage() const;

  // --------------------------- TYPE DATA --------------------------

  /// Special LJ values [0-3] and Special Coul values [4-7]
  UCL_D_Vec<numtyp> sp_lj;

  /// If atom type constants fit in shared memory, use fast kernels
  bool shared_types;

  /// Number of atom types
  int _lj_types;

  numtyp _qqrd2e;

 private:
  bool _allocated;
  numtyp _e_shift, _f_shift, _alpha, _cut_coulsq;
  void loop(const bool _eflag, const bool _vflag);
};

}

#endif
