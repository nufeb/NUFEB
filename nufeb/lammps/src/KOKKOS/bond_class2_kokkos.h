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

#ifdef BOND_CLASS

BondStyle(class2/kk,BondClass2Kokkos<LMPDeviceType>)
BondStyle(class2/kk/device,BondClass2Kokkos<LMPDeviceType>)
BondStyle(class2/kk/host,BondClass2Kokkos<LMPHostType>)

#else

#ifndef LMP_BOND_CLASS2_KOKKOS_H
#define LMP_BOND_CLASS2_KOKKOS_H

#include "bond_class2.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagBondClass2Compute{};

template<class DeviceType>
class BondClass2Kokkos : public BondClass2 {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  BondClass2Kokkos(class LAMMPS *);
  virtual ~BondClass2Kokkos();
  void compute(int, int);
  void coeff(int, char **);
  void read_restart(FILE *);

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondClass2Compute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondClass2Compute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &ebond, const F_FLOAT &fbond, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_x_array_randomread x;
  typename Kokkos::View<double*[3],typename AT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > f;
  typename AT::t_int_2d bondlist;

  Kokkos::DualView<E_FLOAT*,Kokkos::LayoutRight,DeviceType> k_eatom;
  Kokkos::DualView<F_FLOAT*[6],Kokkos::LayoutRight,DeviceType> k_vatom;
  Kokkos::View<E_FLOAT*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_eatom;
  Kokkos::View<F_FLOAT*[6],Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  typename AT::t_ffloat_1d d_k2, d_k3, d_k4;
  typename AT::t_ffloat_1d d_r0;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
