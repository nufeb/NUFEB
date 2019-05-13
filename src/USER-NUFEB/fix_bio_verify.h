/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(verify,FixVerify)

#else

#ifndef LMP_FIX_VERIFY_H
#define LMP_FIX_VERIFY_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixVerify : public Fix {
 public:
  FixVerify(class LAMMPS *, int, char **);
 ~FixVerify();
  void init();
  int setmask();
  void end_of_step();
  int modify_param(int, char **);

 private:

  int nlocal;
  int nall;
  int nnus;                         // # of nutrients

  int nevery;
  int bm1flag, bm2flag, bm3flag, mflag;
  int demflag;

  double **nuS;                    // nutrient concentration for all grids
  double **catCoeff;                 // catabolism coefficients of species
  double **anabCoeff;                // anabolism  coefficients of species
  double **gYield;                   // yield coefficients
  double vol;

  double global_no2, global_pre_no2;
  double global_nh3, global_pre_nh3;
  double global_smass, global_pre_smass;

  class FixKinetics *kinetics;
  class BIO *bio;
  class FixKineticsDiffusion *diffusion;
  class ComputeNufebHeight *cheight;
  class AtomVecBio *avec;

  std::vector< std::vector<int> > nlist;
  int *visit;
  double cutoff;
  std::vector<int> fslist;

  void nitrogen_mass_balance();
  void benchmark_one();
  double get_ave_s_sub_base();
  double get_ave_s_o2_base();
  double get_ave_s_nh4_base();
  void bm1_output();
  void bm3_output();
  void benchmark_two();
  void benchmark_three();
  void remove_atom(double);
  int reachHeight(double);

  // BM2
  void neighbor_list ();
  void free_particle_list();

};

}

#endif
#endif

