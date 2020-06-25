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

FixStyle(kinetics/growth/monod,FixKineticsMonod)

#else

#ifndef SRC_FIX_KINETICSMONOD_H
#define SRC_FIX_KINETICSMONOD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsMonod : public Fix {
 public:
  FixKineticsMonod(class LAMMPS *, int, char **);
  ~FixKineticsMonod();
  void init();
  int setmask();
  void grow_subgrid(int);
  void growth(double, int);

  int external_gflag;

 private:
  double *mu;
  double *decay;
  double *maintain;
  double *yield;
  double **ks;

  int ntypes;

  double *radius;
  double *rmass;
  double *outer_mass;
  double *outer_radius;

  double **nus;
  double **nur;

  double **xdensity;

  int isub, io2, inh4, ino2, ino3, isuc, ico2, ico2g;  // nutrient index
  int ieps;			    // eps index

  int *species;                     // species index 0 = unknown, 1 = het, 2 = aob, 3 = nob, 4 = eps, 5 = dead, 6 = cyano, 7 = ecw
  double ***growrate;               // growth rate [type][epsflag][grid]

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // computational domain size
  int nx, ny, nz;
  double vol;                       // grid volume and gas volume
  double eps_dens;                  // EPS density
  double eta_het;                   // HET reduction factor in anoxic condition
  double suc_exp;                   // Sucrose export rate (0->1)

  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class BIO *bio;

  void init_param();
  void update_biomass(double***, double);

  void growth_het(int, int);
  void growth_aob(int, int);
  void growth_nob(int, int);
  void growth_eps(int, int);
  void growth_dead(int, int);
  void growth_ana(int, int);
  void growth_com(int, int);
  void growth_cyano(int, int);
  void growth_ecw(int,int);
};

}

#endif
#endif
