/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ENABLE_DUMP_BIO_HDF5
#ifdef DUMP_CLASS
DumpStyle(bio/hdf5, DumpBioHDF5)
#else
#ifndef LMP_DUMP_BIO_HDF5_H
#define LMP_DUMP_BIO_HDF5_H

#include "dump.h"

extern "C" {
#include <hdf5.h>
}

#include <string>
#include <vector>

namespace LAMMPS_NS {
// forward declaration
class BIO;
class FixKinetics;
class FixKineticsEnergy;

class DumpBioHDF5 : public Dump {
 public:
  DumpBioHDF5(LAMMPS *, int, char **);
  ~DumpBioHDF5();

protected:
  void write();
  void init_style();
  void write_header(bigint) {}
  void pack(tagint *) {}
  void write_data(int, double *) {}
  int parse_fields(int narg, char **arg);

  hid_t create_filespace_atom(bool oneperproc);
  hid_t create_filespace_grid(bool oneperproc);

  template <class T> 
  herr_t write_atoms_scalar(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset);
  template <class T>
  herr_t write_atoms_comp(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset, int comp);
  template <class T>
  herr_t write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc);
  template <class T>
  herr_t write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int n, int index);

  BIO *bio;
  FixKinetics *kinetics;
  FixKineticsEnergy *energy;
  std::vector<std::string> fields;
};
}

#endif // LMP_DUMP_HDF5_H
#endif // DUMP_CLASS
#endif // ENABLE_DUMP_HDF5
