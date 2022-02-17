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

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_KOKKOS_H
#define LMP_PAIR_KOKKOS_H

#include "Kokkos_Macros.hpp"
#include "pair.h"
#include "neighbor_kokkos.h"
#include "neigh_list_kokkos.h"
#include "Kokkos_Vectorization.hpp"
#include "Kokkos_ScatterView.hpp"

namespace LAMMPS_NS {

template<int Table>
struct CoulLongTable {
  enum {DoTable = Table};
};

// Tags for doing coulomb calculations or not
// They facilitate function overloading, since
// partial template specialization of member functions is not allowed
struct CoulTag {};
struct NoCoulTag {};

template<int FLAG>
struct DoCoul {
  typedef NoCoulTag type;
};

template<>
struct DoCoul<1> {
  typedef CoulTag type;
};


//Specialisation for Neighborlist types Half, HalfThread, Full
template <class PairStyle, int NEIGHFLAG, bool STACKPARAMS, class Specialisation = void>
struct PairComputeFunctor  {
  typedef typename PairStyle::device_type device_type ;
  typedef ArrayTypes<device_type> AT;

  // Reduction type, contains evdwl, ecoul and virial[6]
  typedef EV_FLOAT value_type;

  // The copy of the pair style
  PairStyle c;
  typename AT::t_f_array f;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  // The force array is atomic for Half/Thread neighbor style
  //Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,
  //             device_type,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > f;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,device_type,Kokkos::Experimental::ScatterSum,NeedDup<NEIGHFLAG,device_type>::value > dup_f;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  //Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,
  //             device_type,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > eatom;
  Kokkos::Experimental::ScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,device_type,Kokkos::Experimental::ScatterSum,NeedDup<NEIGHFLAG,device_type>::value > dup_eatom;

  //Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,
  //             device_type,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > vatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,device_type,Kokkos::Experimental::ScatterSum,NeedDup<NEIGHFLAG,device_type>::value > dup_vatom;



  NeighListKokkos<device_type> list;

  PairComputeFunctor(PairStyle* c_ptr,
                          NeighListKokkos<device_type>* list_ptr):
  c(*c_ptr),list(*list_ptr) {
    // allocate duplicated memory
    f = c.f;
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, NeedDup<NEIGHFLAG,device_type>::value >(c.f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, NeedDup<NEIGHFLAG,device_type>::value >(c.d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, NeedDup<NEIGHFLAG,device_type>::value >(c.d_vatom);
  };

  // Call cleanup_copy which sets allocations NULL which are destructed by the PairStyle
  ~PairComputeFunctor() {c.cleanup_copy();list.copymode = 1;};

  KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
    return j >> SBBITS & 3;
  }

  void contribute() {
    Kokkos::Experimental::contribute(c.f, dup_f);

    if (c.eflag_atom)
      Kokkos::Experimental::contribute(c.d_eatom, dup_eatom);

    if (c.vflag_atom)
      Kokkos::Experimental::contribute(c.d_vatom, dup_vatom);
  }

  // Loop over neighbors of one atom without coulomb interaction
  // This function is called in parallel
  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list, const NoCoulTag&) const {

    auto a_f = dup_f.template access<AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;
    const int i = list.d_ilist[ii];
    const X_FLOAT xtmp = c.x(i,0);
    const X_FLOAT ytmp = c.x(i,1);
    const X_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    F_FLOAT fxtmp = 0.0;
    F_FLOAT fytmp = 0.0;
    F_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = neighbors_i(jj);
      const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
      j &= NEIGHMASK;
      const X_FLOAT delx = xtmp - c.x(j,0);
      const X_FLOAT dely = ytmp - c.x(j,1);
      const X_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < c.nlocal)) {
          a_f(j,0) -= delx*fpair;
          a_f(j,1) -= dely*fpair;
          a_f(j,2) -= delz*fpair;
        }

        if (EVFLAG) {
          F_FLOAT evdwl = 0.0;
          if (c.eflag) {
            evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?1.0:0.5)*evdwl;
          }

          if (c.vflag_either || c.eflag_atom) ev_tally(ev,i,j,evdwl,fpair,delx,dely,delz);
        }
      }

    }

    a_f(i,0) += fxtmp;
    a_f(i,1) += fytmp;
    a_f(i,2) += fztmp;

    return ev;
  }

  // Loop over neighbors of one atom with coulomb interaction
  // This function is called in parallel
  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list, const CoulTag& ) const {

    auto a_f = dup_f.template access<AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;
    const int i = list.d_ilist[ii];
    const X_FLOAT xtmp = c.x(i,0);
    const X_FLOAT ytmp = c.x(i,1);
    const X_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);
    const F_FLOAT qtmp = c.q(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    F_FLOAT fxtmp = 0.0;
    F_FLOAT fytmp = 0.0;
    F_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = neighbors_i(jj);
      const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
      const F_FLOAT factor_coul = c.special_coul[sbmask(j)];
      j &= NEIGHMASK;
      const X_FLOAT delx = xtmp - c.x(j,0);
      const X_FLOAT dely = ytmp - c.x(j,1);
      const X_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        F_FLOAT fpair = F_FLOAT();

        if(rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype)))
          fpair+=factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
        if(rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype)))
          fpair+=c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < c.nlocal)) {
          a_f(j,0) -= delx*fpair;
          a_f(j,1) -= dely*fpair;
          a_f(j,2) -= delz*fpair;
        }

        if (EVFLAG) {
          F_FLOAT evdwl = 0.0;
          F_FLOAT ecoul = 0.0;
          if (c.eflag) {
            if(rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype))) {
              evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?1.0:0.5)*evdwl;
            }
            if(rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype))) {
              ecoul = c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);
              ev.ecoul += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?1.0:0.5)*ecoul;
            }
          }

          if (c.vflag_either || c.eflag_atom) ev_tally(ev,i,j,evdwl+ecoul,fpair,delx,dely,delz);
        }
      }
    }

    a_f(i,0) += fxtmp;
    a_f(i,1) += fytmp;
    a_f(i,2) += fztmp;

    return ev;
  }

  // Use TeamPolicy, assume Newton off, Full Neighborlist, and no energy/virial
  // Loop over neighbors of one atom without coulomb interaction
  // This function is called in parallel
  KOKKOS_FUNCTION
  void compute_item_team(Kokkos::TeamPolicy<>::member_type team,
                         const NeighListKokkos<device_type> &list, const NoCoulTag&) const {

    const int inum = team.league_size();
    const int atoms_per_team = team.team_size();
    const int firstatom = team.league_rank()*atoms_per_team;
    const int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const X_FLOAT xtmp = c.x(i,0);
      const X_FLOAT ytmp = c.x(i,1);
      const X_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      t_scalar3<double> fsum;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, t_scalar3<double>& ftmp) {

        int j = neighbors_i(jj);
        const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
        j &= NEIGHMASK;
        const X_FLOAT delx = xtmp - c.x(j,0);
        const X_FLOAT dely = ytmp - c.x(j,1);
        const X_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

          ftmp.x += delx*fpair;
          ftmp.y += dely*fpair;
          ftmp.z += delz*fpair;
        }

      },fsum);

      Kokkos::single(Kokkos::PerThread(team), [&] (){
        f(i,0) += fsum.x;
        f(i,1) += fsum.y;
        f(i,2) += fsum.z;
      });

    });
  }

  // Use TeamPolicy, assume Newton off, Full Neighborlist, and no energy/virial
  // Loop over neighbors of one atom with coulomb interaction
  // This function is called in parallel
  KOKKOS_FUNCTION
  void compute_item_team(Kokkos::TeamPolicy<>::member_type team,
                         const NeighListKokkos<device_type> &list, const CoulTag& ) const {

    const int inum = team.league_size();
    const int atoms_per_team = team.team_size();
    int firstatom = team.league_rank()*atoms_per_team;
    int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const X_FLOAT xtmp = c.x(i,0);
      const X_FLOAT ytmp = c.x(i,1);
      const X_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);
      const F_FLOAT qtmp = c.q(i);

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      t_scalar3<double> fsum;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, t_scalar3<double>& ftmp) {
        int j = neighbors_i(jj);
        const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
        const F_FLOAT factor_coul = c.special_coul[sbmask(j)];
        j &= NEIGHMASK;
        const X_FLOAT delx = xtmp - c.x(j,0);
        const X_FLOAT dely = ytmp - c.x(j,1);
        const X_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          F_FLOAT fpair = F_FLOAT();

          if(rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype)))
            fpair+=factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          if(rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype)))
            fpair+=c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);

          ftmp.x += delx*fpair;
          ftmp.y += dely*fpair;
          ftmp.z += delz*fpair;
        }
      },fsum);

      Kokkos::single(Kokkos::PerThread(team), [&] (){
      f(i,0) += fsum.x;
      f(i,1) += fsum.y;
      f(i,2) += fsum.z;
      });
    });
  }


  // Use TeamPolicy, assume Newton off, Full Neighborlist, and energy/virial
  // Loop over neighbors of one atom without coulomb interaction
  // This function is called in parallel
  KOKKOS_FUNCTION
  EV_FLOAT compute_item_team_ev(Kokkos::TeamPolicy<>::member_type team,
                                const NeighListKokkos<device_type> &list, const NoCoulTag&) const {

    EV_FLOAT ev;

    const int inum = team.league_size();
    const int atoms_per_team = team.team_size();
    const int firstatom = team.league_rank()*atoms_per_team;
    const int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const X_FLOAT xtmp = c.x(i,0);
      const X_FLOAT ytmp = c.x(i,1);
      const X_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      FEV_FLOAT fev;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, FEV_FLOAT& fev_tmp) {

        int j = neighbors_i(jj);
        const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
        j &= NEIGHMASK;
        const X_FLOAT delx = xtmp - c.x(j,0);
        const X_FLOAT dely = ytmp - c.x(j,1);
        const X_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

          fev_tmp.f[0] += delx*fpair;
          fev_tmp.f[1] += dely*fpair;
          fev_tmp.f[2] += delz*fpair;

          F_FLOAT evdwl = 0.0;
          if (c.eflag) {
            evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            fev.evdwl += 0.5*evdwl;
          }
          if (c.vflag_either) {
            fev.v[0] += 0.5*delx*delx*fpair;
            fev.v[1] += 0.5*dely*dely*fpair;
            fev.v[2] += 0.5*delz*delz*fpair;
            fev.v[3] += 0.5*delx*dely*fpair;
            fev.v[4] += 0.5*delx*delz*fpair;
            fev.v[5] += 0.5*dely*delz*fpair;
          }
        }
      },fev);

      Kokkos::single(Kokkos::PerThread(team), [&] (){
        f(i,0) += fev.f[0];
        f(i,1) += fev.f[1];
        f(i,2) += fev.f[2];

        if (c.eflag_global)
          ev.evdwl += fev.evdwl;

        if (c.eflag_atom)
          d_eatom(i,0) += fev.evdwl;

        if (c.vflag_global) {
          ev.v[0] += fev.v[0];
          ev.v[1] += fev.v[1];
          ev.v[2] += fev.v[2];
          ev.v[3] += fev.v[3];
          ev.v[4] += fev.v[4];
          ev.v[5] += fev.v[5];
        }

        if (c.vflag_atom) {
          d_vatom(i,0) += fev.v[0];
          d_vatom(i,1) += fev.v[1];
          d_vatom(i,2) += fev.v[2];
          d_vatom(i,3) += fev.v[3];
          d_vatom(i,4) += fev.v[4];
          d_vatom(i,5) += fev.v[5];
        }
      });
    });
    return ev;
  }

  // Use TeamPolicy, assume Newton off, Full Neighborlist, and energy/virial
  // Loop over neighbors of one atom with coulomb interaction
  // This function is called in parallel
  KOKKOS_FUNCTION
  EV_FLOAT compute_item_team_ev(Kokkos::TeamPolicy<>::member_type team,
                                const NeighListKokkos<device_type> &list, const CoulTag& ) const {

    EV_FLOAT ev;

    const int inum = team.league_size();
    const int atoms_per_team = team.team_size();
    int firstatom = team.league_rank()*atoms_per_team;
    int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const X_FLOAT xtmp = c.x(i,0);
      const X_FLOAT ytmp = c.x(i,1);
      const X_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);
      const F_FLOAT qtmp = c.q(i);

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      FEV_FLOAT fev;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, FEV_FLOAT& fev_tmp) {
        int j = neighbors_i(jj);
        const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
        const F_FLOAT factor_coul = c.special_coul[sbmask(j)];
        j &= NEIGHMASK;
        const X_FLOAT delx = xtmp - c.x(j,0);
        const X_FLOAT dely = ytmp - c.x(j,1);
        const X_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          F_FLOAT fpair = F_FLOAT();

          if(rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype)))
            fpair+=factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          if(rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype)))
            fpair+=c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);

          fev.f[0] += delx*fpair;
          fev.f[1] += dely*fpair;
          fev.f[2] += delz*fpair;

          F_FLOAT evdwl = 0.0;
          F_FLOAT ecoul = 0.0;
          if (c.eflag) {
            if(rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype))) {
              evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              ev.evdwl += 0.5*evdwl;
            }
            if(rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype))) {
              ecoul = c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);
              ev.ecoul += 0.5*ecoul;
            }
          }
          if (c.vflag) {
            fev.v[0] += 0.5*delx*delx*fpair;
            fev.v[1] += 0.5*dely*dely*fpair;
            fev.v[2] += 0.5*delz*delz*fpair;
            fev.v[3] += 0.5*delx*dely*fpair;
            fev.v[4] += 0.5*delx*delz*fpair;
            fev.v[5] += 0.5*dely*delz*fpair;
          }
        }
      },fev);

      Kokkos::single(Kokkos::PerThread(team), [&] (){
        f(i,0) += fev.f[0];
        f(i,1) += fev.f[1];
        f(i,2) += fev.f[2];

        if (c.eflag_global) {
          ev.evdwl += fev.evdwl;
          ev.ecoul += fev.ecoul;
        }

        if (c.eflag_atom)
          d_eatom(i,0) += fev.evdwl + fev.ecoul;

        if (c.vflag_global) {
          ev.v[0] += fev.v[0];
          ev.v[1] += fev.v[1];
          ev.v[2] += fev.v[2];
          ev.v[3] += fev.v[3];
          ev.v[4] += fev.v[4];
          ev.v[5] += fev.v[5];
        }

        if (c.vflag_atom) {
          d_vatom(i,0) += fev.v[0];
          d_vatom(i,1) += fev.v[1];
          d_vatom(i,2) += fev.v[2];
          d_vatom(i,3) += fev.v[3];
          d_vatom(i,4) += fev.v[4];
          d_vatom(i,5) += fev.v[5];
        }
      });
    });
    return ev;
  }

  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const
  {
    auto a_eatom = dup_eatom.template access<AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_vatom = dup_vatom.template access<AtomicDup<NEIGHFLAG,device_type>::value>();

    const int EFLAG = c.eflag;
    const int NEWTON_PAIR = c.newton_pair;
    const int VFLAG = c.vflag_either;

    if (EFLAG) {
      if (c.eflag_atom) {
        const E_FLOAT epairhalf = 0.5 * epair;
        if (NEWTON_PAIR || i < c.nlocal) a_eatom[i] += epairhalf;
        if ((NEWTON_PAIR || j < c.nlocal) && NEIGHFLAG != FULL) a_eatom[j] += epairhalf;
      }
    }

    if (VFLAG) {
      const E_FLOAT v0 = delx*delx*fpair;
      const E_FLOAT v1 = dely*dely*fpair;
      const E_FLOAT v2 = delz*delz*fpair;
      const E_FLOAT v3 = delx*dely*fpair;
      const E_FLOAT v4 = delx*delz*fpair;
      const E_FLOAT v5 = dely*delz*fpair;

      if (c.vflag_global) {
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR) {
            ev.v[0] += v0;
            ev.v[1] += v1;
            ev.v[2] += v2;
            ev.v[3] += v3;
            ev.v[4] += v4;
            ev.v[5] += v5;
          } else {
            if (i < c.nlocal) {
              ev.v[0] += 0.5*v0;
              ev.v[1] += 0.5*v1;
              ev.v[2] += 0.5*v2;
              ev.v[3] += 0.5*v3;
              ev.v[4] += 0.5*v4;
              ev.v[5] += 0.5*v5;
            }
            if (j < c.nlocal) {
              ev.v[0] += 0.5*v0;
              ev.v[1] += 0.5*v1;
              ev.v[2] += 0.5*v2;
              ev.v[3] += 0.5*v3;
              ev.v[4] += 0.5*v4;
              ev.v[5] += 0.5*v5;
            }
          }
        } else {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
      }

      if (c.vflag_atom) {
        if (NEWTON_PAIR || i < c.nlocal) {
          a_vatom(i,0) += 0.5*v0;
          a_vatom(i,1) += 0.5*v1;
          a_vatom(i,2) += 0.5*v2;
          a_vatom(i,3) += 0.5*v3;
          a_vatom(i,4) += 0.5*v4;
          a_vatom(i,5) += 0.5*v5;
        }
        if ((NEWTON_PAIR || j < c.nlocal) && NEIGHFLAG != FULL) {
          a_vatom(j,0) += 0.5*v0;
          a_vatom(j,1) += 0.5*v1;
          a_vatom(j,2) += 0.5*v2;
          a_vatom(j,3) += 0.5*v3;
          a_vatom(j,4) += 0.5*v4;
          a_vatom(j,5) += 0.5*v5;
        }
      }
    }
  }


  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (c.newton_pair) compute_item<0,1>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
    else compute_item<0,0>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &energy_virial) const {
    if (c.newton_pair)
      energy_virial += compute_item<1,1>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
    else
      energy_virial += compute_item<1,0>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<>::member_type& team) const {
    compute_item_team(team,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<>::member_type& team, value_type &energy_virial) const {
    energy_virial += compute_item_team_ev(team,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }
};

template <class PairStyle, bool STACKPARAMS, class Specialisation>
struct PairComputeFunctor<PairStyle,N2,STACKPARAMS,Specialisation>  {
  typedef typename PairStyle::device_type device_type ;
  typedef EV_FLOAT value_type;

  PairStyle c;
  NeighListKokkos<device_type> list;

  PairComputeFunctor(PairStyle* c_ptr,
                          NeighListKokkos<device_type>* list_ptr):
  c(*c_ptr),list(*list_ptr) {};
  ~PairComputeFunctor() {c.cleanup_copy();list.copymode = 1;};

  KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
    return j >> SBBITS & 3;
  }


  void contribute() {}

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list, const NoCoulTag&) const {
    (void) list;
    EV_FLOAT ev;
    const int i = ii;//list.d_ilist[ii];
    const X_FLOAT xtmp = c.x(i,0);
    const X_FLOAT ytmp = c.x(i,1);
    const X_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);

    //const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = c.nall;

    F_FLOAT fxtmp = 0.0;
    F_FLOAT fytmp = 0.0;
    F_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jj;//neighbors_i(jj);
      if(i==j) continue;
      const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
      j &= NEIGHMASK;
      const X_FLOAT delx = xtmp - c.x(j,0);
      const X_FLOAT dely = ytmp - c.x(j,1);
      const X_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        if (EVFLAG) {
          F_FLOAT evdwl = 0.0;
          if (c.eflag) {
            evdwl = 0.5*
              factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            ev.evdwl += evdwl;
          }

          if (c.vflag_either || c.eflag_atom) ev_tally(ev,i,j,evdwl,fpair,delx,dely,delz);
        }
      }
    }

    c.f(i,0) += fxtmp;
    c.f(i,1) += fytmp;
    c.f(i,2) += fztmp;

    return ev;
  }

  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const
  {
    const int EFLAG = c.eflag;
    const int VFLAG = c.vflag_either;

    if (EFLAG) {
      if (c.eflag_atom) {
        const E_FLOAT epairhalf = 0.5 * epair;
        if (i < c.nlocal) c.d_eatom[i] += epairhalf;
        if (j < c.nlocal) c.d_eatom[j] += epairhalf;
      }
    }

    if (VFLAG) {
      const E_FLOAT v0 = delx*delx*fpair;
      const E_FLOAT v1 = dely*dely*fpair;
      const E_FLOAT v2 = delz*delz*fpair;
      const E_FLOAT v3 = delx*dely*fpair;
      const E_FLOAT v4 = delx*delz*fpair;
      const E_FLOAT v5 = dely*delz*fpair;

      if (c.vflag_global) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
      }

      if (c.vflag_atom) {
        if (i < c.nlocal) {
          c.d_vatom(i,0) += 0.5*v0;
          c.d_vatom(i,1) += 0.5*v1;
          c.d_vatom(i,2) += 0.5*v2;
          c.d_vatom(i,3) += 0.5*v3;
          c.d_vatom(i,4) += 0.5*v4;
          c.d_vatom(i,5) += 0.5*v5;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    compute_item<0,0>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &energy_virial) const {
    energy_virial += compute_item<1,0>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<>::member_type& team) const
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<>::member_type& team, value_type &energy_virial) const
  {}

};

// Filter out Neighflags which are not supported for PairStyle
// The enable_if clause will invalidate the last parameter of the function, so that
// a match is only achieved, if PairStyle supports the specific neighborlist variant.
// This uses the fact that failure to match template parameters is not an error.
// By having the enable_if with a ! and without it, exactly one of the functions
// pair_compute_neighlist will match - either the dummy version
// or the real one further below.
template<class PairStyle, unsigned NEIGHFLAG, class Specialisation>
EV_FLOAT pair_compute_neighlist (PairStyle* fpair, typename Kokkos::Impl::enable_if<!((NEIGHFLAG&PairStyle::EnabledNeighFlags) != 0), NeighListKokkos<typename PairStyle::device_type>*>::type list) {
  EV_FLOAT ev;
  (void) fpair;
  (void) list;
  printf("ERROR: calling pair_compute with invalid neighbor list style: requested %i  available %i \n",NEIGHFLAG,PairStyle::EnabledNeighFlags);
  return ev;
}

template<class FunctorStyle>
int GetTeamSize(FunctorStyle& functor, int team_size, int vector_length) {
    int team_size_max = Kokkos::TeamPolicy<>::team_size_max(functor);

#ifdef KOKKOS_ENABLE_CUDA
    if(team_size*vector_length > team_size_max)
      team_size = team_size_max/vector_length;
#else
    team_size = 1;
#endif
    return team_size;
}

// Submit ParallelFor for NEIGHFLAG=HALF,HALFTHREAD,FULL,N2
template<class PairStyle, unsigned NEIGHFLAG, class Specialisation>
EV_FLOAT pair_compute_neighlist (PairStyle* fpair, typename Kokkos::Impl::enable_if<(NEIGHFLAG&PairStyle::EnabledNeighFlags) != 0, NeighListKokkos<typename PairStyle::device_type>*>::type list) {
  EV_FLOAT ev;

  if (!fpair->lmp->kokkos->neigh_thread_set)
    if (list->inum <= 16384 && NEIGHFLAG == FULL)
      fpair->lmp->kokkos->neigh_thread = 1;

  if (fpair->lmp->kokkos->neigh_thread) {
    int vector_length = 8;
    int atoms_per_team = 32;

    if(fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
      PairComputeFunctor<PairStyle,NEIGHFLAG,false,Specialisation > ff(fpair,list);
      atoms_per_team = GetTeamSize(ff, atoms_per_team, vector_length);
      Kokkos::TeamPolicy<Kokkos::IndexType<int> > policy(list->inum,atoms_per_team,vector_length);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(policy,ff,ev);
      else                              Kokkos::parallel_for(policy,ff);
    } else {
      PairComputeFunctor<PairStyle,NEIGHFLAG,true,Specialisation > ff(fpair,list);
      atoms_per_team = GetTeamSize(ff, atoms_per_team, vector_length);
      Kokkos::TeamPolicy<Kokkos::IndexType<int> > policy(list->inum,atoms_per_team,vector_length);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(policy,ff,ev);
      else                              Kokkos::parallel_for(policy,ff);
    }
  } else {
    if(fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
      PairComputeFunctor<PairStyle,NEIGHFLAG,false,Specialisation > ff(fpair,list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else                              Kokkos::parallel_for(list->inum,ff);
      ff.contribute();
    } else {
      PairComputeFunctor<PairStyle,NEIGHFLAG,true,Specialisation > ff(fpair,list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else                              Kokkos::parallel_for(list->inum,ff);
      ff.contribute();
    }
  }
  return ev;
}

template<class PairStyle, class Specialisation>
EV_FLOAT pair_compute (PairStyle* fpair, NeighListKokkos<typename PairStyle::device_type>* list) {
  EV_FLOAT ev;
  if (fpair->neighflag == FULL) {
    ev = pair_compute_neighlist<PairStyle,FULL,Specialisation> (fpair,list);
  } else if (fpair->neighflag == HALFTHREAD) {
    ev = pair_compute_neighlist<PairStyle,HALFTHREAD,Specialisation> (fpair,list);
  } else if (fpair->neighflag == HALF) {
    ev = pair_compute_neighlist<PairStyle,HALF,Specialisation> (fpair,list);
  } else if (fpair->neighflag == N2) {
    ev = pair_compute_neighlist<PairStyle,N2,Specialisation> (fpair,list);
  }
  return ev;
}

template<class DeviceType>
struct PairVirialFDotRCompute {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typename AT::t_x_array_const_um x;
  typename AT::t_f_array_const_um f;
  const int offset;

  PairVirialFDotRCompute(  typename AT::t_x_array_const_um x_,
  typename AT::t_f_array_const_um f_,
  const int offset_):x(x_),f(f_),offset(offset_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int j, value_type &energy_virial) const {
    const int i = j + offset;
    energy_virial.v[0] += f(i,0)*x(i,0);
    energy_virial.v[1] += f(i,1)*x(i,1);
    energy_virial.v[2] += f(i,2)*x(i,2);
    energy_virial.v[3] += f(i,1)*x(i,0);
    energy_virial.v[4] += f(i,2)*x(i,0);
    energy_virial.v[5] += f(i,2)*x(i,1);
  }
};

template<class PairStyle>
void pair_virial_fdotr_compute(PairStyle* fpair) {
  EV_FLOAT virial;
  if (fpair->neighbor->includegroup == 0) {
    int nall = fpair->atom->nlocal + fpair->atom->nghost;
    Kokkos::parallel_reduce(nall,PairVirialFDotRCompute<typename PairStyle::device_type>(fpair->x,fpair->f,0),virial);
  } else {
    Kokkos::parallel_reduce(fpair->atom->nfirst,PairVirialFDotRCompute<typename PairStyle::device_type>(fpair->x,fpair->f,0),virial);
    EV_FLOAT virial_ghost;
    Kokkos::parallel_reduce(fpair->atom->nghost,PairVirialFDotRCompute<typename PairStyle::device_type>(fpair->x,fpair->f,fpair->atom->nlocal),virial_ghost);
    virial+=virial_ghost;
  }
  fpair->vflag_fdotr = 0;
  fpair->virial[0] = virial.v[0];
  fpair->virial[1] = virial.v[1];
  fpair->virial[2] = virial.v[2];
  fpair->virial[3] = virial.v[3];
  fpair->virial[4] = virial.v[4];
  fpair->virial[5] = virial.v[5];
}




}

#endif
#endif

/* ERROR/WARNING messages:

*/
