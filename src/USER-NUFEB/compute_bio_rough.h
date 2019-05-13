/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(roughness,ComputeNufebRough)

#else

#ifndef LMP_COMPUTE_ROUGH_H
#define LMP_COMPUTE_ROUGH_H

#include "compute.h"
#include "reduce_grid.h"

namespace LAMMPS_NS {

class ComputeNufebRough : public Compute, public ReduceGrid<ComputeNufebRough> {
  friend ReduceGrid<ComputeNufebRough>;

 public:
  ComputeNufebRough(class LAMMPS *, int, char **);
  virtual ~ComputeNufebRough() {}
  void init();
  void grow_subgrid();
  double compute_scalar();

 private:
  int nx, ny;
  double stepx, stepy;
  std::vector<double> maxh;
  Grid<double, 2> grid;
  Subgrid<double, 2> subgrid;

  Subgrid<double, 2> get_subgrid() const { return subgrid; }
  template <typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      *result++ = maxh[*it];
    }
    return result;
  }
  template <typename InputIterator0, typename InputIterator1, typename BinaryOperation>
  InputIterator1 unpack_cells_reduce(InputIterator0 first, InputIterator0 last, InputIterator1 input, BinaryOperation op) {
    for (InputIterator0 it = first; it != last; ++it) {
      maxh[*it] = op(maxh[*it], *input++);
    }
    return input;
  }
  int get_cell_data_size(int n) const { return n; }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
