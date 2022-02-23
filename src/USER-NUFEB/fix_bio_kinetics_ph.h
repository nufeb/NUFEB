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

FixStyle(kinetics/ph,FixKineticsPH)

#else

#ifndef SRC_FIX_KINETICSPH_H
#define SRC_FIX_KINETICSPH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsPH : public Fix {
 public:
  FixKineticsPH(class LAMMPS *, int, char **);
  ~FixKineticsPH();
  void init();
  int setmask();
  void solve_ph();
  void buffer_ph();

  double buffer_flag;              // 1 = buffer ph, 0 = unbuffer ph

 private:
  class FixKinetics *kinetics;
  class BIO *bio;

  double phflag;                   // 0 = fix ph, 1 = dynamic ph
  double **keq;                    // equilibrium constants [nutrient][4]
  double iph;                      // initial ph
  double phlo, phhi;               // lower and upper bounds of ph buffer

  void output_data();
  void compute_activity(int, int, double);
  void init_keq();
  void dynamic_ph(int, int);
};

}

#endif
#endif

