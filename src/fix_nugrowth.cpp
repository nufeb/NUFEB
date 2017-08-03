/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <USER-NUFEB/fix_nugrowth.h>
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixNuGrowth::FixNuGrowth(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 20) error->all(FLERR,"Illegal fix growth command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix growth command");

  var = new char*[16];
  ivar = new int[16];

  int i;
  for (i = 0; i < 16; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

//  pFile[0] = fopen ("Mass HET 60", "w");
//  pFile[1] = fopen ("Mass AOB 60", "w");
//  pFile[2] = fopen ("Mass NOB 60", "w");
}

/* ---------------------------------------------------------------------- */

FixNuGrowth::~FixNuGrowth()
{
  int i;
  for (i = 0; i < 16; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixNuGrowth::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNuGrowth::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix growth requires atom attribute diameter");

  int i;
  for (i = 0; i < 16; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix nugrowth does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix nugrowth is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void FixNuGrowth::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  change_dia();
}

/* ---------------------------------------------------------------------- */
unsigned int int_to_int(unsigned int k) {
    return (k == 0 || k == 1 ? k : ((k % 2) + 10 * int_to_int(k / 2)));
}

void FixNuGrowth::change_dia()
{

  modify->clearstep_compute();

  double KsHET = input->variable->compute_equal(ivar[0]);
  double Ko2HET = input->variable->compute_equal(ivar[1]);
  double Kno2HET = input->variable->compute_equal(ivar[2]);
  double Kno3HET = input->variable->compute_equal(ivar[3]);
  double Knh4AOB = input->variable->compute_equal(ivar[4]);
  double Ko2AOB = input->variable->compute_equal(ivar[5]);
  double Kno2NOB = input->variable->compute_equal(ivar[6]);
  double Ko2NOB = input->variable->compute_equal(ivar[7]);
  double MumHET = input->variable->compute_equal(ivar[8]);
  double MumAOB = input->variable->compute_equal(ivar[9]);
  double MumNOB = input->variable->compute_equal(ivar[10]);
  double etaHET = input->variable->compute_equal(ivar[11]);
  // double bHET = input->variable->compute_equal(ivar[12]);
  // double bAOB = input->variable->compute_equal(ivar[13]);
  // double bNOB = input->variable->compute_equal(ivar[14]);
  double bEPS = input->variable->compute_equal(ivar[12]);
  double YEPS = input->variable->compute_equal(ivar[13]);
  double YHET = input->variable->compute_equal(ivar[14]);
  double EPSdens = input->variable->compute_equal(ivar[15]);

  double density;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outerMass = atom->outerMass;
  double *outerRadius = atom->outerRadius;
  double *sub = atom->sub;
  double *o2 = atom->o2;
  double *nh4 = atom->nh4;
  double *no2 = atom->no2;
  double *no3 = atom->no3;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double R1 = 0.0;
  double R2 = 0.0;
  double R3 = 0.0;
  double R4 = 0.0;
  double R5 = 0.0;

  int i;
  for (i = 0; i < nall; i++) {
//	  fprintf(stdout, "groubit =%i\n",groupbit);
//	  fprintf(stdout, "mask =%i\n",int_to_int(mask[i]));
//	  fprintf(stdout, "type =%i\n",type[i]);

    if (mask[i] & groupbit) {
      double gHET = 0;
      double gAOB = 0;
      double gNOB = 0;
      double gEPS = 0;
      if (type[i] == 1) {
        gHET = 1;
      }
      if (type[i] == 2) {
        gAOB = 1;
      }
      if (type[i] == 3) {
        gNOB = 1;
      }
      if (type[i] == 4) {
        gEPS = 1;
      }

      R1 = MumHET*(sub[i]/(KsHET+sub[i]))*(o2[i]/(Ko2HET+o2[i]));
      R2 = MumAOB*(nh4[i]/(Knh4AOB+nh4[i]))*(o2[i]/(Ko2AOB+o2[i]));
      R3 = MumNOB*(no2[i]/(Kno2NOB+no2[i]))*(o2[i]/(Ko2NOB+o2[i]));
      R4 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no3[i]/(Kno3HET+no3[i]))*(Ko2HET/(Ko2HET+o2[i]));
      R5 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no2[i]/(Kno2HET+no2[i]))*(Ko2HET/(Ko2HET+o2[i]));
      // double R6 = bHET;
      // double R7 = bAOB;
      // double R8 = bNOB;
      double R9 = bEPS;

      double value = update->dt * (gHET*(R1+R4+R5) + gAOB*R2 + gNOB*R3 - gEPS*R9);

      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      double oldMass = rmass[i];
      rmass[i] = rmass[i]*(1 + (value*nevery));
      if (rmass[i] <= 0) {
        rmass[i] = oldMass;
      }
      // fprintf(stdout, "Radius: %e Outer Radius: %e\n", radius[i], outerRadius[i]);
      // fprintf(stdout, "ID: %i Type: %i Outer Mass: %e\n", atom->tag[i], atom->type[i], outerMass[i]);
      

      double value2 = update->dt * (YEPS/YHET)*(R1+R4+R5);
      outerMass[i] = (((4.0*MY_PI/3.0)*((outerRadius[i]*outerRadius[i]*outerRadius[i])-(radius[i]*radius[i]*radius[i])))*EPSdens)+(value2*nevery*rmass[i]);

      outerRadius[i] = pow((3.0/(4.0*MY_PI))*((rmass[i]/density)+(outerMass[i]/EPSdens)),(1.0/3.0));

      radius[i] = pow(((6*rmass[i])/(density*MY_PI)),(1.0/3.0))*0.5;
    }
  }
//
//  //output analytical total mass value
//  if(update->ntimestep%10 == 0){
//		double amhet = 1.0e-16*exp((R1+R4+R5)*update->ntimestep*60);
//		double amaob = 1.0e-16*exp((R2)*update->ntimestep*60);
//		double amnob = 1.0e-16*exp((R3)*update->ntimestep*60);
//
//		double mhet = compute_totalmass(1);
//	    double maob = compute_totalmass(2);
//	    double mnob = compute_totalmass(3);
//
//		int nhet = compute_totaln(1);
//		int naob = compute_totaln(2);
//		int nnob = compute_totaln(3);
//		fprintf(pFile[0], "%i\t %e\t %e\t %i\n",update->ntimestep, mhet, amhet, nhet);
//		fprintf(pFile[1], "%i\t %e\t %e\t %i\n",update->ntimestep, maob, amaob, naob);
//		fprintf(pFile[2], "%i\t %e\t %e\t %i\n",update->ntimestep, mnob, amnob, nnob);
//  }

  modify->addstep_compute(update->ntimestep + nevery);
}


double FixNuGrowth::compute_totalmass(int type){
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double tmass = 0;

  int i;
  for (i = 0; i < nall; i++) {
	if (atom->type[i] == type) {
		tmass += atom->rmass[i];
	}
  }
  //fprintf(stdout, "%e\n",tmass);
  return tmass;
}

int FixNuGrowth::compute_totaln(int type){
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int tnum = 0;

  int i;
  for (i = 0; i < nall; i++) {
	if (atom->type[i] == type) {
		tnum ++;
	}
  }
  //fprintf(stdout, "%e\n",tmass);
  return tnum;
}
