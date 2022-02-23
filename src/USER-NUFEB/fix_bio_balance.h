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

FixStyle(kinetics/balance,FixKineticsBalance)

#else

#ifndef SRC_FIX_KINETICS_BALANCE_H
#define SRC_FIX_KINETICS_BALANCE_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsBalance : public Fix {
 public:
  FixKineticsBalance(class LAMMPS *, int, char **);
  ~FixKineticsBalance();
  int setmask();
  void post_constructor();
  void init();
  void setup(int);
  void setup_pre_exchange();
  void pre_exchange();
  void pre_neighbor();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  int nevery,lbstyle,nitermax;
  double thresh,stopthresh;
  char bstr[4];
  int wtflag;                   // 1 for weighted balancing

  double imbnow;                // current imbalance factor
  double imbprev;               // imbalance factor before last rebalancing
  double imbfinal;              // imbalance factor after last rebalancing
  double maxloadperproc;        // max load on any processor
  int itercount;                // iteration count of last call to Balance
  int kspace_flag;              // 1 if KSpace solver defined
  int pending;
  bigint lastbalance;           // last timestep balancing was attempted

  class FixKinetics *kinetics;
  class Balance *balance;
  class Irregular *irregular;

  void rebalance();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix balance shift string is invalid

The string can only contain the characters "x", "y", or "z".

E: Fix balance rcb cannot be used with comm_style brick

Comm_style tiled must be used instead.

E: Cannot open fix balance output file

Self-explanatory.

*/
