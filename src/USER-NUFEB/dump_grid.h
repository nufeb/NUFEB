/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef ENABLE_DUMP_GRID
#ifdef DUMP_CLASS
DumpStyle(grid, DumpGrid)
#else
#ifndef LMP_DUMP_GRID_H
#define LMP_DUMP_GRID_H

#include "dump.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include <functional>
#include <vector>

namespace LAMMPS_NS {
// forward declaration
class BIO;
class FixKinetics;
class FixKineticsEnergy;

class DumpGrid : public Dump {
 public:
  DumpGrid(LAMMPS *, int, char **);
  ~DumpGrid();

 protected:
  void write();
  void init_style();
  void write_header(bigint) {}
  void pack(tagint *) {}
  void write_data(int, double *) {}
  int parse_fields(int narg, char **arg);

  void pack_concentration(vtkSmartPointer<vtkImageData>);
  void pack_uptake(vtkSmartPointer<vtkImageData>);
  void pack_activity(vtkSmartPointer<vtkImageData>);
  void pack_yield(vtkSmartPointer<vtkImageData>);
  void pack_catabolism(vtkSmartPointer<vtkImageData>);
  void pack_anabolism(vtkSmartPointer<vtkImageData>);
  void pack_hydronium(vtkSmartPointer<vtkImageData>);
  void pack_tuple1(vtkSmartPointer<vtkImageData>, const char *, double *);
  void pack_tuple1(vtkSmartPointer<vtkImageData>, const char *, double **, char **, int);
  void pack_tuple5(vtkSmartPointer<vtkImageData>, const char *, double ***, char **, int);

  BIO *bio;
  FixKinetics *kinetics;
  FixKineticsEnergy *energy;
  std::vector<std::string> fields;
  std::vector<std::function<void(vtkSmartPointer<vtkImageData>)>> packs;
};
}

#endif // LMP_DUMP_GRID_H
#endif // DUMP_CLASS
#endif // ENABLE_DUMP_GRID
