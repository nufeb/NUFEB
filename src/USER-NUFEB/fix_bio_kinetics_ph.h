/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

