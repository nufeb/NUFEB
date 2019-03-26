/*
 * ibm.cpp
 *
 *  Created on: 8 Nov 2016
 *      Author: bowen
 */

#include "bio.h"

#include <string.h>
#include <cctype>
#include <cstdio>
#include <cstdlib>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "atom_vec_bio.h"

using namespace LAMMPS_NS;

BIO::BIO(LAMMPS *lmp) : Pointers(lmp)
{
  //type
  yield = NULL;
  maintain = NULL;
  decay = NULL;
  q = NULL;
  mu = NULL;
  ks = NULL;
  edoner = NULL;
  tname = NULL;
  anab_coeff = NULL;
  cata_coeff = NULL;
  decay_coeff = NULL;
  tgibbs_coeff = NULL;
  dissipation = NULL;
  tgibbs_flag = NULL;
  tcharge = NULL;

  //nutrient
  nnu = 0;
  ini_nus = NULL;
  nuname = NULL;
  diff_coeff = NULL;
  nustate = NULL;
  nugibbs_coeff = NULL;
  ngflag = NULL;
  nucharge = NULL;
  kla = NULL;
  mw = NULL;
}

/* ---------------------------------------------------------------------- */

BIO::~BIO()
{
  memory->destroy(yield);
  memory->destroy(maintain);
  memory->destroy(decay);
  memory->destroy(edoner);
  memory->destroy(q);
  memory->destroy(mu);
  memory->destroy(ks);
  memory->destroy(anab_coeff);
  memory->destroy(cata_coeff);
  memory->destroy(decay_coeff);
  memory->destroy(tgibbs_coeff);
  memory->destroy(dissipation);
  memory->destroy(tgibbs_flag);
  memory->destroy(tcharge);

  memory->destroy(diff_coeff);
  memory->destroy(mw);
  memory->destroy(ini_nus);
  memory->destroy(nustate);
  memory->destroy(nugibbs_coeff);
  memory->destroy(ngflag);
  memory->destroy(nucharge);
  memory->destroy(kla);

  for (int i = 0; i < atom->ntypes+1; i++) {
    delete [] tname[i];
  }

  for (int i = 0; i < nnu+1; i++) {
    delete [] nuname[i];
  }

  memory->sfree(tname);
  memory->sfree(nuname);
}

void BIO::type_grow()
{
  if (yield != NULL) memory->grow(yield,atom->ntypes+1,"bio:yield");
  if (maintain != NULL) memory->grow(maintain,atom->ntypes+1,"bio:maintain");
  if (decay != NULL) memory->grow(decay,atom->ntypes+1,"bio:decay");
  if (dissipation != NULL) memory->grow(dissipation,atom->ntypes+1,"bio:dissipation");
  if (q != NULL) memory->grow(q,atom->ntypes+1,"bio:q");
  if (mu != NULL) memory->grow(mu,atom->ntypes+1,"bio:mu");
  if (ks != NULL) memory->grow(ks,atom->ntypes+1,nnu+1,"bio:ks");
  if (anab_coeff != NULL) memory->grow(anab_coeff,atom->ntypes+1,nnu+1,"bio:anab_coeff");
  if (cata_coeff != NULL) memory->grow(cata_coeff,atom->ntypes+1,nnu+1,"bio:cata_coeff");
  if (decay_coeff != NULL) memory->grow(decay_coeff,atom->ntypes+1,nnu+1,"bio:decay_coeff");
  if (tgibbs_coeff != NULL) memory->grow(tgibbs_coeff,atom->ntypes+1,5,"bio:tgibbs_coeff");
  if (tgibbs_flag != NULL) memory->grow(tgibbs_flag,atom->ntypes+1,"bio:tgflag");
  if (tcharge != NULL) memory->grow(tcharge,atom->ntypes+1,5,"bio:tcharge");
  if (edoner != NULL) memory->grow(edoner,atom->ntypes+1,"bio:edoner");
  if (tname != NULL) tname = (char **) memory->srealloc(tname,(atom->ntypes+1)*sizeof(char *),"bio:tname");
}

void BIO::create_type(char *name) {
  atom->ntypes = atom->ntypes + 1;
  type_grow();
  tname[atom->ntypes] = name;
}

void BIO::data_nutrients(int narg, char **arg)
{
  if (narg != 10) error->all(FLERR,"Incorrect args for nutrient definitions");

  int id = force->numeric(FLERR,arg[0]);
  double scell = force->numeric(FLERR,arg[3]);
  double xbcm = force->numeric(FLERR,arg[4]);
  double xbcp = force->numeric(FLERR,arg[5]);
  double ybcm = force->numeric(FLERR,arg[6]);
  double ybcp = force->numeric(FLERR,arg[7]);
  double zbcm = force->numeric(FLERR,arg[8]);
  double zbcp = force->numeric(FLERR,arg[9]);

  //nutrient name
  char *name;
  int n = strlen(arg[1]) + 1;
  name = new char[n];
  strcpy(name,arg[1]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(name[i]) && name[i] != '_')
      error->all(FLERR,"Nutrient name must be "
                 "alphanumeric or underscore characters");

  for (int i = 0; i < nnu+1; i++)
    if ((nuname[i] != NULL) && (strcmp(nuname[i], name) == 0)
        && (i != id)){
      error->one(FLERR,"Repeat nutrient names");
    }

  if (nuname[id] == NULL) {
    nuname[id] = new char[n];
  } else if (strcmp(nuname[id], name) != 0){
    error->one(FLERR,"Incompatible nutrient names");
  }

  strcpy(nuname[id],name);

  delete [] name;

  //nutrient type
  int m = strlen(arg[2]);
  if (m != 1) error->all(FLERR,"Nutrient type must be a single char, "
      "l = liq, g = gas");
  char type = arg[2][0];
  if (type == 'l') nustate[id] = 0;
  else if (type == 'g') nustate[id] = 1;
  else error->all(FLERR,"Undefined nutrient type, "
      "l = liq, g = gas");

  if (ini_nus == NULL) error->all(FLERR,"Cannot set nutrient concentration for this nutrient style");
  ini_nus[id][0] = scell;
  ini_nus[id][1] = xbcm;
  ini_nus[id][2] = xbcp;
  ini_nus[id][3] = ybcm;
  ini_nus[id][4] = ybcp;
  ini_nus[id][5] = zbcm;
  ini_nus[id][6] = zbcp;
}

void BIO::set_tname(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for type name definitions");
  int id = force->numeric(FLERR,arg[0]);
  //type name
  char *name;
  int n = strlen(arg[1]) + 1;
  name = new char[n];
  strcpy(name,arg[1]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(name[i]) && name[i] != '_')
      error->all(FLERR,"Type name must be "
                       "alphanumeric or underscore characters");

  for (int i = 0; i < atom->ntypes; i++)
    if ((tname[i] != NULL) && (strcmp(tname[i], name) == 0)
        && (i != id)){
      error->one(FLERR,"Duplicated type name");
    }

  if (tname[id] == NULL) {
    tname[id] = new char[n];
  } else if (strcmp(tname[id], name) != 0){
    error->one(FLERR,"Incompatible type name");
  }

  strcpy(tname[id],name);

  delete [] name;
}

/* ----------------------------------------------------------------------
   set uptake rate for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_q(const char *str)
{
  if (q == NULL) error->all(FLERR,"Cannot set consumption rate for this atom style");

  char* name;
  double growth_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&growth_one);

  if (n != 2) error->all(FLERR,"Invalid consumption rate line in data file");

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for consumption rate set");

  q[itype] = growth_one;

  if (q[itype] < 0.0) error->all(FLERR,"Invalid consumption rate value");
}

/* ----------------------------------------------------------------------
   set growth values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_mu(const char *str)
{
  if (mu == NULL) error->all(FLERR,"Cannot set growth rate for this atom style");

  char* name;
  double growth_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&growth_one);

  if (n != 2) error->all(FLERR,"Invalid growth line in data file");

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for growth set");

  mu[itype] = growth_one;

  if (mu[itype] < 0.0) error->all(FLERR,"Invalid growth value");

  AtomVecBio* avec = (AtomVecBio *) atom->style_match("bio");
}

/* ----------------------------------------------------------------------
   set ks values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_ks(int narg, char **arg)
{
  if (ks == NULL) error->all(FLERR,"Cannot set Ks for this atom style");
  if (narg != nnu+1) error->all(FLERR,"Invalid Ks line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for Ks set");

  for(int i = 1; i < nnu+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    ks[itype][i] = value;
  }
}

/* ----------------------------------------------------------------------
   set yield values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_yield(const char *str)
{
  if (yield == NULL) error->all(FLERR,"Cannot set yield for this atom style");

  char* name;
  double yield_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&yield_one);
  if (n != 2) error->all(FLERR,"Invalid set_yield line in data file");

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_yield set");

  if (yield_one < 0)
    lmp->error->all(FLERR,"yield cannot be zero or less than zero");

  yield[itype] = yield_one;

  if (yield[itype] < 0.0) error->all(FLERR,"Invalid set_yield value");
}


/* ----------------------------------------------------------------------
   set maintenance values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_maintain(const char *str)
{
  if (maintain == NULL) error->all(FLERR,"Cannot set maintain for this atom style");

  char* name;
  double maintain_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&maintain_one);
  if (n != 2) error->all(FLERR,"Invalid set_maintain line in data file");

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_maintain set");

  maintain[itype] = maintain_one;
  //mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set decay values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_decay(const char *str)
{
  if (decay == NULL) error->all(FLERR,"Cannot set decay for this atom style");

  char* name;
  double decay_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&decay_one);
  if (n != 2) error->all(FLERR,"Invalid set_decay line in data file");

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_decay set");

  decay[itype] = decay_one;
  //mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set diffusion coeff for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_diffusion(const char *str)
{
  if (diff_coeff == NULL) error->all(FLERR,"Cannot set diffCoeff for this nutrient style");

  char* name;
  double diffu_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&diffu_one);
  if (n != 2) error->all(FLERR,"Invalid diffCoeff line in data file");

  int inu = find_nuid(name);
  delete [] name;

  if (inu < 1 || inu > nnu)
    error->all(FLERR,"Invalid nutrient for diffusion coefficient set");

  diff_coeff[inu] = diffu_one;
  //mass_setflag[itype] = 1;

  if (diff_coeff[inu] < 0.0) error->all(FLERR,"Invalid diffCoeff value");
}

/* ----------------------------------------------------------------------
   set molecular weights for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_mw(const char *str)
{
  if (mw == NULL) error->all(FLERR,"Cannot set molecular weights for this nutrient style");

  char* name;
  double mw_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&mw_one);
  if (n != 2) error->all(FLERR,"Invalid molecular weights line in data file");

  int inu = find_nuid(name);
  delete [] name;

  if (inu < 1 || inu > nnu)
    error->all(FLERR,"Invalid nutrient for molecular weights set");

  mw[inu] = mw_one;

  if (mw[inu] < 0.0) error->all(FLERR,"Invalid molecular weights value");
}


/* ----------------------------------------------------------------------
   set dissipation values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_dissipation(const char *str)
{
  if (dissipation == NULL) error->all(FLERR,"Cannot set dissipation for this atom style");

  char* name;
  double diss_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&diss_one);
  if (n != 2) error->all(FLERR,"Invalid dissipation line in data file");

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for dissipation set");

  dissipation[itype] = diss_one;

  if (dissipation[itype] < 0.0) error->all(FLERR,"Invalid dissipation value");
}

/* ----------------------------------------------------------------------
   set catabolism coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_cata_coeff(int narg, char **arg)
{
  if (cata_coeff == NULL) error->all(FLERR,"Cannot set catCoeff for this atom style");
  if (narg != nnu+1) error->all(FLERR,"Invalid catCoeff line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for catabolism coefficient set");

  for(int i = 1; i < nnu+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    cata_coeff[itype][i] = value;
  }
}

/* ----------------------------------------------------------------------
   set anabolism coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_anab_coeff(int narg, char **arg)
{
  if (anab_coeff == NULL) error->all(FLERR,"Cannot set anabCoeff for this atom style");
  if (narg != nnu+1) error->all(FLERR,"Invalid anabCoeff line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for anabolism coefficient set");

  for(int i = 1; i < nnu+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    anab_coeff[itype][i] = value;
  }
}

/* ----------------------------------------------------------------------
   set decay coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_decay_coeff(int narg, char **arg)
{
  if (decay_coeff == NULL) error->all(FLERR,"Cannot set decayCoeff for this atom style");
  if (narg != nnu+1) error->all(FLERR,"Invalid decayCoeff line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for decay coefficient set");

  for(int i = 1; i < nnu+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    decay_coeff[itype][i] = value;
  }
}

/* ----------------------------------------------------------------------
   set energy coefficient for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_nugibbs_coeff(int narg, char **arg)
{
  if (nugibbs_coeff == NULL) error->all(FLERR,"Cannot set energy coeff for this nutrient");
  if (narg != 7) error->all(FLERR,"Invalid nuGCOeff line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int inu = find_nuid(name);
  delete [] name;

  if (inu < 1 || inu > nnu)
    error->all(FLERR,"Invalid nutrient for nuG coefficient set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "inf") == 0) {
      nugibbs_coeff[inu][i] = INF;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      nugibbs_coeff[inu][i] = value;
    }
  }

  int flag = atoi(arg[6]);

  if ((flag > 0) && (flag < 6) && (nugibbs_coeff[inu][flag-1] != INF))
    ngflag[inu] = flag - 1;
  else error->all(FLERR,"Invalid nutrient energy flag");
}

/* ----------------------------------------------------------------------
   set energy coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_tgibbs_coeff(int narg, char **arg)
{
  if (tgibbs_coeff == NULL) error->all(FLERR,"Cannot set energy coeff for this type");
  if (narg != 7) error->all(FLERR,"Invalid typeGCoeff line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for typeG coefficient set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "inf") == 0) {
      tgibbs_coeff[itype][i] = INF;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      tgibbs_coeff[itype][i] = value;
    }
  }

  int flag = atoi(arg[6]);
  if ((flag > 0) && (flag < 6) && (tgibbs_coeff[itype][flag-1] != INF))
    tgibbs_flag[itype] = flag - 1;
  else error->all(FLERR,"Invalid type energy flag");
}

/* ----------------------------------------------------------------------
   set Electron Donor values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_edoner(int narg, char **arg)
{
  if (edoner == NULL) error->all(FLERR,"Cannot set Electron Donor for this atom style");
  if (narg != 2) error->all(FLERR,"Incorrect Electron Donor for type name definitions");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  edoner[itype] = -1;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for Electron Donor set");

  for(int i = 1; i < nnu+1; i++) {
    if (strcmp(arg[1], nuname[i]) == 0) {
      edoner[itype] = i;
    } else if (strcmp(arg[1], "null") == 0) {
      edoner[itype] = 0;
    }
  }

  if (edoner[itype] < 0) error->all(FLERR,"Cannot find Electron Donor for this atom style");
}

/* ----------------------------------------------------------------------
   set charges for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_nucharge(int narg, char **arg)
{
  if (nucharge == NULL) error->all(FLERR,"Cannot set charge for this nutrient");
  if (narg != 6) error->all(FLERR,"Invalid nutrient charge line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int inu = find_nuid(name);
  delete [] name;

  if (inu < 1 || inu > nnu)
    error->all(FLERR,"Invalid nutrient for nutrient charge set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "na") == 0) {
      nucharge[inu][i] = NA;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      nucharge[inu][i] = value;
    }
  }
}

/* ----------------------------------------------------------------------
   set charges for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_tcharge(int narg, char **arg)
{
  if (tcharge == NULL) error->all(FLERR,"Cannot set charge for this type");
  if (narg != 6) error->all(FLERR,"Invalid type charge line in data file");

  char* name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  int itype = find_typeid(name);
  delete [] name;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for typeG coefficient set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "na") == 0) {
      tcharge[itype][i] = 10001;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      tcharge[itype][i] = value;
    }
  }
}

/* ----------------------------------------------------------------------
   set mass transfer coeff for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_kla(const char *str)
{
  if (kla == NULL) error->all(FLERR,"Cannot set KLa for this atom style");

  char* name;
  double kLa_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str,"%s %lg",name,&kLa_one);
  if (n != 2) error->all(FLERR,"Invalid KLa line in data file");

  int inu = find_nuid(name);
  delete [] name;

  if (inu < 1 || inu > nnu)
    error->all(FLERR,"Invalid nutrient for KLa set");

  kla[inu] = kLa_one;
  //mass_setflag[itype] = 1;

  if (kla[inu] < 0.0) error->all(FLERR,"Invalid KLa value");
}

int BIO::find_typeid(char *name) {

  for (int i = 1; i < atom->ntypes+1; i++){
    if (tname[i] && strcmp(tname[i],name) == 0) {
      return i;
    }
  }
  return -1;
}

int BIO::find_nuid(char *name) {

  for (int i = 1; i < nnu+1; i++)
    if (nuname[i] && strcmp(nuname[i],name) == 0) {
      return i;
    }

  return -1;
}
