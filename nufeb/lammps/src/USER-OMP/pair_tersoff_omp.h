/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(tersoff/omp,PairTersoffOMP)

#else

#ifndef LMP_PAIR_TERSOFF_OMP_H
#define LMP_PAIR_TERSOFF_OMP_H

#include "pair_tersoff.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairTersoffOMP : public PairTersoff, public ThrOMP {

 public:
  PairTersoffOMP(class LAMMPS *);

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  template <int EVFLAG, int EFLAG, int VFLAG_ATOM>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif
