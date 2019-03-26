/*
 * bio.h
 *
 *  Created on: 8 Nov 2016
 *      Author: bowen
 */

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

  double *diff_coeff;         // diffusion coefficient [nutrient]
  double *mw;                 // molecular Weights [nutrient]
  double **ini_nus;           // inlet nutrient concentrations [nutrient][1grid + 5bc]
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
