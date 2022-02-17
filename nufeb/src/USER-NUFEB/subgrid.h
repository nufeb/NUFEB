/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_SUBGRID_H
#define LMP_SUBGRID_H

#include "grid.h"

#include <cmath>
#include <functional>
#include <vector>

namespace LAMMPS_NS {
namespace detail {
template <typename T, std::size_t N, std::size_t I>
struct for_ {
  static void loop(const std::array<T, N> &dimensions,
		   std::array<T, N> &index,
		   std::function<void(const std::array<T, N> &)> f) {
    for (int i = 0; i < dimensions[I - 1]; i++) {
      index[I - 1] = i;
      for_<T, N, I - 1>::loop(dimensions, index, f); 
    }
  }
};
 
template <typename T, std::size_t N>
struct for_<T, N, 0> {
  static void loop(const std::array<T, N> &dimensions,
		   std::array<T, N> &index,
		   std::function<void(const std::array<T, N> &)> f) {
    f(index);
  }
};
} // namespace detail

template <typename T, std::size_t N>
struct for_ {
  static void loop(const std::array<T, N> &dimensions,
		   std::function<void(const std::array<T, N> &)> f) {
    std::array<T, N> index;
    detail::for_<T, N, N>::loop(dimensions, index, f);
  }
};
 
template <typename T, std::size_t N, typename IndexType = int>
class Subgrid {
 public:
  Subgrid() = default;
  Subgrid(const Subgrid &other) = default;
  Subgrid(const Grid<T, N, IndexType> &grid, const std::array<IndexType, N> &origin, const std::array<IndexType, N> &dimensions) :
    grid(grid), origin(origin), dimensions(dimensions) {}
  Subgrid(const Grid<T, N, IndexType> &grid, const Box<IndexType, N> & box) :
    Subgrid(grid, box.lower, size(box)) {}
  Subgrid(const Grid<T, N, IndexType> &grid,
	  const Box<T, N> & box,
	  std::function<T (T)> round = [](T value) { return std::floor(value); }) :
    grid(grid) {
    for (int i = 0; i < N; i++) {
      origin[i] = static_cast<IndexType>(round(box.lower[i] / grid.get_cell_size()[i]));
      dimensions[i] = static_cast<IndexType>(round(box.upper[i] / grid.get_cell_size()[i])) - origin[i];
    }
  }

  const Grid<T, N, IndexType> &get_grid() const {
    return grid;
  }

  const std::array<IndexType, N> &get_origin() const {
    return origin;
  }
 
  const std::array<IndexType, N> &get_dimensions() const {
    return dimensions;
  }

  Box<IndexType, N> get_box() const {
    Box<IndexType, N> result;
    result.lower = origin;
    for (int i = 0; i < N; i++) {
      result.upper[i] = origin[i] + dimensions[i];
    }
    return result;
  }

  IndexType get_linear_index(const std::array<IndexType, N> &cell) const {
    IndexType n = T(1);
    IndexType result = T(0);
    for (int i = 0; i < N; i++) {
      result += (cell[i] - origin[i]) * n;
      n *= dimensions[i];
    }
    return result;
  }

  IndexType cell_count() const {
    IndexType result = T(1);
    for(int i = 0; i < N; i++) {
      result *= dimensions[i];
    }
    return result;
  }

  IndexType get_index(const std::array<T, N> &x) const {
    IndexType n = T(1);
    IndexType result = T(0);
    for (int i = 0; i < N; i++) {
      result += (static_cast<IndexType>(x[i] / grid.get_cell_size()[i]) - origin[i]) * n;
      n *= dimensions[i];
    }
    return result;
  }

  bool is_inside(const std::array<T, N> &x) const {
    bool result = true;
    for (int i = 0; i < N; i++) {
      IndexType index = static_cast<IndexType>(x[i] / grid.get_cell_size()[i]);
      result &= index >= origin[i] && index < origin[i] + dimensions[i];  
    }
    return result;
  }

  // TODO: enable this only if T is_floating_point
  // note that this will be subject to RVO (copy elision)
  std::vector<std::array<T, N>> get_cell_centers() const {
    std::vector<std::array<T, N>> result;
    for_<IndexType, N>::loop(dimensions,
	    [&](const std::array<IndexType, N> &index) {
	      std::array<T, N> x;
	      for (int i = 0; i < N; i++) {
		x[i] = grid.get_cell_size()[i] / 2 + (origin[i] + index[i]) * grid.get_cell_size()[i];
	      }
	      result.push_back(x);
	    });
    return result;
  }

 private:
  Grid<T, N, IndexType> grid;
  std::array<IndexType, N> origin;
  std::array<IndexType, N> dimensions;
};
} // namespace LAMMPS_NS
#endif // LMP_SUBGRID_H
