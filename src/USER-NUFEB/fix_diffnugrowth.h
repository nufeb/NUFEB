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

#ifdef FIX_CLASS

FixStyle(diffnugrowth,FixDiffNuGrowth)

#else

#ifndef LMP_FIX_DIFFNUGROWTH_H
#define LMP_FIX_DIFFNUGROWTH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDiffNuGrowth : public Fix {
 public:
  FixDiffNuGrowth(class LAMMPS *, int, char **);
  ~FixDiffNuGrowth();
  int setmask();
  void init();
  void pre_force(int);

 private:

  char **var;
  int *ivar;
  int diffevery;
  int outputevery;
  double *xCell;
  double *yCell;
  double *zCell;
  double *cellVol;
  bool *ghost;
  double *subCell;
  double *o2Cell;
  double *nh4Cell;
  double *no2Cell;
  double *no3Cell;

  double KsHET;
  double Ko2HET;
  double Kno2HET;
  double Kno3HET;
  double Knh4AOB;
  double Ko2AOB;
  double Kno2NOB;
  double Ko2NOB;
  double MumHET;
  double MumAOB;
  double MumNOB;
  double etaHET;
  double bHET; // R6
  double bAOB; // R7
  double bNOB; // R8
  double bEPS; // R9
  double bmHET;
  double bmAOB;
  double bmNOB;
  double bX;
  double YHET;
  double YAOB;
  double YNOB;
  double YEPS;
  double Y1;
  double EPSdens;
  double Do2;
  double Dnh4;
  double Dno2;
  double Dno3;
  double Ds;
  double diffT;

  int numCells;
  int nx, ny, nz;
  double initsub, inito2, initnh4, initno2, initno3;
  double subBC, o2BC, no2BC, no3BC, nh4BC;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int bflag; // 1 = dirichlet, 2 = neumann, 3 = mixed
  double xstep, ystep, zstep;
  double sumRs, sumRo2, sumRno2, sumRno3, sumRnh4;
  void change_dia();
  void compute_flux(double *, double *, double *, double, double, double, int);
  void output_data(int, int);
  bool is_convergence(double *, double *, double, double);
  int overlap();
  void compute_Rvalues(double*, double*, double*, double*, double*);
};

}

#endif
#endif
