/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pair_snap_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "force.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "kokkos.h"
#include "sna.h"

#define MAXLINE 1024
#define MAXWORD 3

namespace LAMMPS_NS {

// Outstanding issues with quadratic term
// 1. there seems to a problem with compute_optimized energy calc
// it does not match compute_regular, even when quadratic coeffs = 0

//static double t1 = 0.0;
//static double t2 = 0.0;
//static double t3 = 0.0;
//static double t4 = 0.0;
//static double t5 = 0.0;
//static double t6 = 0.0;
//static double t7 = 0.0;
/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkos<DeviceType>::PairSNAPKokkos(LAMMPS *lmp) : PairSNAP(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_cutsq = tdual_fparams("PairSNAPKokkos::cutsq",atom->ntypes+1,atom->ntypes+1);
  auto d_cutsq = k_cutsq.template view<DeviceType>();
  rnd_cutsq = d_cutsq;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkos<DeviceType>::~PairSNAPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->request(this,instance_me);

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == HALF || neighflag == HALFTHREAD) { // still need atomics, even though using a full neigh list
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else {
    error->all(FLERR,"Must use half neighbor list style with pair snap/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& max_neighs) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }
};

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  copymode = 1;
  int newton_pair = force->newton_pair;
  if (newton_pair == false)
    error->all(FLERR,"PairSNAPKokkos requires 'newton on'");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  /*
  for (int i = 0; i < nlocal; i++) {
    typename t_neigh_list::t_neighs neighs_i = neigh_list.get_neighs(i);
    const int num_neighs = neighs_i.get_num_neighs();
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }*/
  max_neighs = 0;
  Kokkos::parallel_reduce("PairSNAPKokkos::find_max_neighs",inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Experimental::Max<int>(max_neighs));

  int vector_length = 1;
  int team_size = 1;
  int team_size_max = Kokkos::TeamPolicy<DeviceType>::team_size_max(*this);
#ifdef KOKKOS_ENABLE_CUDA
  team_size = 32;//max_neighs;
  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
#endif

  if (beta_max < inum) {
    beta_max = inum;
    d_beta = Kokkos::View<F_FLOAT**, DeviceType>("PairSNAPKokkos:beta",inum,ncoeff);
    d_ninside = Kokkos::View<int*, DeviceType>("PairSNAPKokkos:ninside",inum);
  }

  int chunk_size = MIN(2000,inum);
  chunk_offset = 0;

  snaKK.grow_rij(chunk_size,max_neighs);

  EV_FLOAT ev;

  while (chunk_offset < inum) { // chunk up loop to prevent running out of memory

    EV_FLOAT ev_tmp;

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    //ComputeNeigh
    typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeNeigh> policy_neigh(chunk_size,team_size,vector_length);
    Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);

    //PreUi
    typename Kokkos::RangePolicy<DeviceType, TagPairSNAPPreUi> policy_preui(0,chunk_size);
    Kokkos::parallel_for("PreUi",policy_preui,*this);

    //ComputeUi
    typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUi> policy_ui(((inum+team_size-1)/team_size)*max_neighs,team_size,vector_length);
    Kokkos::parallel_for("ComputeUi",policy_ui,*this);

    //Ulisttot transpose
    snaKK.transpose_ulisttot();

    //Compute bispectrum
    if (quadraticflag || eflag) {
      //ComputeZi
      int idxz_max = snaKK.idxz_max;
      typename Kokkos::RangePolicy<DeviceType, TagPairSNAPComputeZi> policy_zi(0,chunk_size*idxz_max);
      Kokkos::parallel_for("ComputeZi",policy_zi,*this);

      //ComputeBi
      typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeBi> policy_bi(chunk_size,team_size,vector_length);
      Kokkos::parallel_for("ComputeBi",policy_bi,*this);
    }

    //Compute beta = dE_i/dB_i for all i in list
    typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPBeta> policy_beta(chunk_size,team_size,vector_length);
    Kokkos::parallel_for("ComputeBeta",policy_beta,*this);

    //ComputeYi
    typename Kokkos::RangePolicy<DeviceType, TagPairSNAPZeroYi> policy_zero_yi(0,chunk_size);
    Kokkos::parallel_for("ZeroYi",policy_zero_yi,*this);

    //ComputeYi
    int idxz_max = snaKK.idxz_max;
    typename Kokkos::RangePolicy<DeviceType, TagPairSNAPComputeYi> policy_yi(0,chunk_size*idxz_max);
    Kokkos::parallel_for("ComputeYi",policy_yi,*this);

    //ComputeDuidrj
    typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDuidrj> policy_duidrj(((inum+team_size-1)/team_size)*max_neighs,team_size,vector_length);
    Kokkos::parallel_for("ComputeDuidrj",policy_duidrj,*this);

    //ComputeDeidrj
    typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDeidrj> policy_deidrj(((inum+team_size-1)/team_size)*max_neighs,team_size,vector_length);
    Kokkos::parallel_for("ComputeDeidrj",policy_deidrj,*this);

    //ComputeForce
    if (eflag) {
      if (neighflag == HALF) {
        typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<HALF,1> > policy_force(chunk_size,team_size,vector_length);
        Kokkos::parallel_reduce(policy_force
          ,*this,ev_tmp);
      } else if (neighflag == HALFTHREAD) {
        typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<HALFTHREAD,1> > policy_force(chunk_size,team_size,vector_length);
        Kokkos::parallel_reduce(policy_force
          ,*this,ev_tmp);
      }
    } else {
      if (neighflag == HALF) {
        typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<HALF,0> > policy_force(chunk_size,team_size,vector_length);
        Kokkos::parallel_for(policy_force
          ,*this);
      } else if (neighflag == HALFTHREAD) {
        typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<HALFTHREAD,0> > policy_force(chunk_size,team_size,vector_length);
        Kokkos::parallel_for(policy_force
          ,*this);
      }
    }
    ev += ev_tmp;
    chunk_offset += chunk_size;
  } // end while

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  atomKK->modified(execution_space,F_MASK);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f     = decltype(dup_f)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPBeta,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPBeta>::member_type& team) const {

  // TODO: use RangePolicy instead, or thread over ncoeff?
  int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const int ielem = d_map[itype];
  SNAKokkos<DeviceType> my_sna = snaKK;

  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    d_coeffi(d_coeffelem,ielem,Kokkos::ALL);

  for (int icoeff = 0; icoeff < ncoeff; icoeff++)
    d_beta(ii,icoeff) = d_coeffi[icoeff+1];

  if (quadraticflag) {
    int k = ncoeff+1;
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      double bveci = my_sna.blist(ii,icoeff);
      d_beta(ii,icoeff) += d_coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
        double bvecj = my_sna.blist(ii,jcoeff);
        d_beta(ii,icoeff) += d_coeffi[k]*bvecj;
        d_beta(ii,jcoeff) += d_coeffi[k]*bveci;
        k++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::allocate()
{
  PairSNAP::allocate();

  int n = atom->ntypes;
  d_map = Kokkos::View<T_INT*, DeviceType>("PairSNAPKokkos::map",n+1);
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairSNAPKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairSNAP::init_one(i,j);
  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairSNAP::coeff(narg,arg);

  // Set up element lists

  d_radelem = Kokkos::View<F_FLOAT*, DeviceType>("pair:radelem",nelements);
  d_wjelem = Kokkos::View<F_FLOAT*, DeviceType>("pair:wjelem",nelements);
  d_coeffelem = Kokkos::View<F_FLOAT**, Kokkos::LayoutRight, DeviceType>("pair:coeffelem",nelements,ncoeffall);

  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_wjelem = Kokkos::create_mirror_view(d_wjelem);
  auto h_coeffelem = Kokkos::create_mirror_view(d_coeffelem);
  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int ielem = 0; ielem < nelements; ielem++) {
    h_radelem(ielem) = radelem[ielem];
    h_wjelem(ielem) = wjelem[ielem];
    for (int jcoeff = 0; jcoeff < ncoeffall; jcoeff++) {
      h_coeffelem(ielem,jcoeff) = coeffelem[ielem][jcoeff];
    }
  }

  for (int i = 1; i <= atom->ntypes; i++) {
    h_map(i) = map[i];
  }

  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_wjelem,h_wjelem);
  Kokkos::deep_copy(d_coeffelem,h_coeffelem);
  Kokkos::deep_copy(d_map,h_map);

  snaKK = SNAKokkos<DeviceType>(rfac0,twojmax,
                  rmin0,switchflag,bzeroflag);
  snaKK.grow_rij(0,0);
  snaKK.init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeNeigh>::member_type& team) const {

  int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  SNAKokkos<DeviceType> my_sna = snaKK;
  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);
  const int itype = type[i];
  const int ielem = d_map[itype];
  const double radi = d_radelem[ielem];

  const int num_neighs = d_numneigh[i];

  // rij[][3] = displacements between atom I and those neighbors
  // inside = indices of neighbors of I within cutoff
  // wj = weights for neighbors of I within cutoff
  // rcutij = cutoffs for neighbors of I within cutoff
  // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

  int ninside = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,num_neighs),
      [&] (const int jj, int& count) {
    Kokkos::single(Kokkos::PerThread(team), [&] (){
      T_INT j = d_neighbors(i,jj);
      const F_FLOAT dx = x(j,0) - xtmp;
      const F_FLOAT dy = x(j,1) - ytmp;
      const F_FLOAT dz = x(j,2) - ztmp;

      const int jtype = type(j);
      const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
      const int elem_j = d_map[jtype];

      if ( rsq < rnd_cutsq(itype,jtype) )
       count++;
    });
  },ninside);

  d_ninside(ii) = ninside;

  if (team.team_rank() == 0)
  Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_neighs),
      [&] (const int jj, int& offset, bool final) {
  //for (int jj = 0; jj < num_neighs; jj++) {
    T_INT j = d_neighbors(i,jj);
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;

    const int jtype = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
    const int elem_j = d_map[jtype];

    if ( rsq < rnd_cutsq(itype,jtype) ) {
      if (final) {
        my_sna.rij(ii,offset,0) = dx;
        my_sna.rij(ii,offset,1) = dy;
        my_sna.rij(ii,offset,2) = dz;
        my_sna.inside(ii,offset) = j;
        my_sna.wj(ii,offset) = d_wjelem[elem_j];
        my_sna.rcutij(ii,offset) = (radi + d_radelem[elem_j])*rcutfac;
      }
      offset++;
    }
  });
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPPreUi,const int& ii) const {
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.pre_ui(ii);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeUi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeUi>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((inum+team.team_size()-1)/team.team_size()));
  if (ii >= inum) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((inum+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_ui(team,ii,jj);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPZeroYi,const int& ii) const {
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.zero_yi(ii);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeYi,const int& ii) const {
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.compute_yi(ii,d_beta);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeZi,const int& ii) const {
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.compute_zi(ii);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeBi,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeBi>::member_type& team) const {
  int ii = team.league_rank();
  SNAKokkos<DeviceType> my_sna = snaKK;
  my_sna.compute_bi(team,ii);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeDuidrj,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDuidrj>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((inum+team.team_size()-1)/team.team_size()));
  if (ii >= inum) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((inum+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_duidrj(team,ii,jj);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeDeidrj,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeDeidrj>::member_type& team) const {
  SNAKokkos<DeviceType> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((inum+team.team_size()-1)/team.team_size()));
  if (ii >= inum) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((inum+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_deidrj(team,ii,jj);
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  SNAKokkos<DeviceType> my_sna = snaKK;
  const int ninside = d_ninside(ii);

  Kokkos::parallel_for (Kokkos::TeamThreadRange(team,ninside),
      [&] (const int jj) {
    int j = my_sna.inside(ii,jj);

    F_FLOAT fij[3];
    fij[0] = my_sna.dedr(ii,jj,0);
    fij[1] = my_sna.dedr(ii,jj,1);
    fij[2] = my_sna.dedr(ii,jj,2);

    Kokkos::single(Kokkos::PerThread(team), [&] (){
      a_f(i,0) += fij[0];
      a_f(i,1) += fij[1];
      a_f(i,2) += fij[2];
      a_f(j,0) -= fij[0];
      a_f(j,1) -= fij[1];
      a_f(j,2) -= fij[2];

      // tally global and per-atom virial contribution

      if (EVFLAG) {
        if (vflag_either) {
          v_tally_xyz<NEIGHFLAG>(ev,i,j,
            fij[0],fij[1],fij[2],
            -my_sna.rij(ii,jj,0),-my_sna.rij(ii,jj,1),
            -my_sna.rij(ii,jj,2));
        }
      }

    });
  });

  // tally energy contribution

  if (EVFLAG) {
    if (eflag_either) {

      const int itype = type(i);
      const int ielem = d_map[itype];
      Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        d_coeffi(d_coeffelem,ielem,Kokkos::ALL);

      Kokkos::single(Kokkos::PerTeam(team), [&] () {

        // evdwl = energy of atom I, sum over coeffs_k * Bi_k

        double evdwl = d_coeffi[0];

        // E = beta.B + 0.5*B^t.alpha.B

        // linear contributions

        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          evdwl += d_coeffi[icoeff+1]*my_sna.blist(ii,icoeff);

        // quadratic contributions

        if (quadraticflag) {
          int k = ncoeff+1;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bveci = my_sna.blist(ii,icoeff);
            evdwl += 0.5*d_coeffi[k++]*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              double bvecj = my_sna.blist(ii,jcoeff);
              evdwl += d_coeffi[k++]*bveci*bvecj;
            }
          }
        }

        //ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
        if (eflag_global) ev.evdwl += evdwl;
        if (eflag_atom) d_eatom[i] += evdwl;
      });
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG> >::member_type& team) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>(), team, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz,
      const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const E_FLOAT v0 = delx*fx;
  const E_FLOAT v1 = dely*fy;
  const E_FLOAT v2 = delz*fz;
  const E_FLOAT v3 = delx*fy;
  const E_FLOAT v4 = delx*fz;
  const E_FLOAT v5 = dely*fz;

  if (vflag_global) {
    ev.v[0] += v0;
    ev.v[1] += v1;
    ev.v[2] += v2;
    ev.v[3] += v3;
    ev.v[4] += v4;
    ev.v[5] += v5;
  }

  if (vflag_atom) {
    a_vatom(i,0) += 0.5*v0;
    a_vatom(i,1) += 0.5*v1;
    a_vatom(i,2) += 0.5*v2;
    a_vatom(i,3) += 0.5*v3;
    a_vatom(i,4) += 0.5*v4;
    a_vatom(i,5) += 0.5*v5;
    a_vatom(j,0) += 0.5*v0;
    a_vatom(j,1) += 0.5*v1;
    a_vatom(j,2) += 0.5*v2;
    a_vatom(j,3) += 0.5*v3;
    a_vatom(j,4) += 0.5*v4;
    a_vatom(j,5) += 0.5*v5;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

template<class DeviceType>
double PairSNAPKokkos<DeviceType>::memory_usage()
{
  double bytes = Pair::memory_usage();
  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);
  bytes += n*n*sizeof(double);
  bytes += (2*ncoeffall)*sizeof(double);
  bytes += (ncoeff*3)*sizeof(double);
  bytes += snaKK.memory_usage();
  return bytes;
}

}
