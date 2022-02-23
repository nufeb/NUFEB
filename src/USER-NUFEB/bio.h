/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef SRC_IBM_H_
#define SRC_IBM_H_

#include "pointers.h"

namespace LAMMPS_NS {

enum {
  INF = INT_MAX,
  NA = 0,
};

class BIO : protected Pointers {
 public:
  //type (microbial species)
  char **tname;               // type name
  double **ks;                // half-saturation constant [type][nutrient]
  double *q;                  // specific consumption rate
  double *mu;                 // specific growth rate
  double *yield;              // growth yield coefficient
  double *dissipation;        // universal gas constant (thermodynamics)
  double *maintain;           // maintenance [type]
  double *decay;              // decay rate [type]
  int *edoner;                // electron donor [type]

  double **cata_coeff;        // catabolism coefficient [type][nutrient]
  double **anab_coeff;        // anabolism coefficient [type][nutrient]
  double **decay_coeff;       // decay coefficient [type][nutrient]
  double **tgibbs_coeff;      // Gibbs free energy coefficient [type][5charges]
  int *tgibbs_flag;           // Gibbs free energy flag for type [type]
  int **tcharge;              // charge [type][5charges]

  //nutrient
  int nnu;                    // # of nutrients
  int *nustate;               // nutrient types 0 = liq, 1 = gas
  char **nuname;              // nutrient name
  int **nubc;                 // boundary condition type [nutrient][3surface pairs]

  double *diff_coeff;         // diffusion coefficient [nutrient]
  double *mw;                 // molecular Weights [nutrient]
  double **init_nus;           // inlet nutrient concentrations [nutrient][1grid + 1bc]
  double **nugibbs_coeff;     // Gibbs free energy coefficient [nutrient][5charges]
  int *ngflag;                // Gibbs free energy flag for nutrients [nutrient]
  int **nucharge;             // charge [nutrient][5charges]
  double *kla;                // mass Transfer Coefficient [nutrient]

  BIO(class LAMMPS *);
  ~BIO();

  void type_grow();
  void create_type(char *);
  void data_nutrients(int, char **);
  void set_tname(int narg, char **arg);
  void set_q(const char *);
  void set_mu(const char *);
  void set_mw(const char *);
  void set_ks(int, char **);
  void set_yield(const char *);
  void set_edoner(int narg, char **arg);
  void set_maintain(const char *);
  void set_decay(const char *);
  void set_diffusion(const char *);
  void set_cata_coeff(int, char **);
  void set_decay_coeff(int, char **);
  void set_anab_coeff(int, char **);
  void set_nugibbs_coeff(int, char **);
  void set_tgibbs_coeff(int, char **);
  void set_dissipation(const char *);
  void set_nucharge(int, char **);
  void set_tcharge(int, char **);
  void set_kla(const char *);
  void set_group_mask();

  int find_typeid(char *name);
  int find_nuid(char *name);

};

}

#endif /* SRC_IBM_H_ */
