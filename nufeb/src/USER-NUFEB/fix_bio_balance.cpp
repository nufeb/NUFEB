/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_bio_balance.h"
#include "fix_bio_kinetics.h"
#include "balance.h"
#include "update.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "irregular.h"
#include "force.h"
#include "kspace.h"
#include "modify.h"
#include "fix_store.h"
#include "rcb.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SHIFT,BISECTION};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

FixKineticsBalance::FixKineticsBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), balance(NULL), irregular(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix balance command");

  box_change_domain = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
  global_freq = 1;

  // parse required arguments

  int dimension = domain->dimension;

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix balance command");
  thresh = force->numeric(FLERR,arg[4]);

  if (strcmp(arg[5],"shift") == 0) lbstyle = SHIFT;
  else if (strcmp(arg[5],"rcb") == 0) lbstyle = BISECTION;
  else error->all(FLERR,"Illegal fix balance command");

  int iarg = 5;
  if (lbstyle == SHIFT) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal fix balance command");
    if (strlen(arg[iarg+1]) > 3)
      error->all(FLERR,"Illegal fix balance command");
    strcpy(bstr,arg[iarg+1]);
    nitermax = force->inumeric(FLERR,arg[iarg+2]);
    if (nitermax <= 0) error->all(FLERR,"Illegal fix balance command");
    stopthresh = force->numeric(FLERR,arg[iarg+3]);
    if (stopthresh < 1.0) error->all(FLERR,"Illegal fix balance command");
    iarg += 4;
  } else if (lbstyle == BISECTION) {
    iarg++;
  }

  // error checks

  if (lbstyle == SHIFT) {
    int blen = strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Fix balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Fix balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Fix balance shift string is invalid");
    }
  }

  if (lbstyle == BISECTION && comm->style == 0)
    error->all(FLERR,"Fix balance rcb cannot be used with comm_style brick");

  // create instance of Balance class
  // if SHIFT, initialize it with params
  // process remaining optional args via Balance

  balance = new Balance(lmp);
  if (lbstyle == SHIFT) balance->shift_setup(bstr,nitermax,thresh);
  balance->options(iarg,narg,arg);
  wtflag = balance->wtflag;

  if (balance->varflag && nevery == 0)
    error->all(FLERR,"Fix balance nevery = 0 cannot be used with weight var");

  // create instance of Irregular class

  irregular = new Irregular(lmp);

  // only force reneighboring if nevery > 0

  if (nevery) force_reneighbor = 1;
  lastbalance = -1;
  
  // compute initial outputs

  itercount = 0;
  pending = 0;
  imbfinal = imbprev = maxloadperproc = 0.0;
}

/* ---------------------------------------------------------------------- */

FixKineticsBalance::~FixKineticsBalance()
{
  delete balance;
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixKineticsBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsBalance::post_constructor()
{
  if (wtflag) balance->weight_storage(id);
}

/* ---------------------------------------------------------------------- */

void FixKineticsBalance::init()
{
  kinetics = NULL;
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  balance->init_imbalance(1);
}

/* ---------------------------------------------------------------------- */

void FixKineticsBalance::setup(int vflag)
{
  // compute final imbalance factor if setup_pre_exchange() invoked balancer
  // this is called at end of run setup, before output

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixKineticsBalance::setup_pre_exchange()
{
  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // insure atoms are in current box & update box via shrink-wrap
  // has to be be done before rebalance() invokes Irregular::migrate_atoms()
  //   since it requires atoms be inside simulation box
  //   even though pbc() will be done again in Verlet::run()
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshhold exceeded

  balance->set_weights();
  imbnow = balance->imbalance_factor(maxloadperproc);
  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixKineticsBalance::pre_exchange()
{
  // return if not a rebalance timestep

  if (nevery && update->ntimestep < next_reneighbor) return;

  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshhold exceeded
  // if weight variable is used, wrap weight setting in clear/add compute

  if (balance->varflag) modify->clearstep_compute();
  balance->set_weights();
  if (balance->varflag) modify->addstep_compute(update->ntimestep + nevery);

  imbnow = balance->imbalance_factor(maxloadperproc);
  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   compute final imbalance factor based on nlocal after comm->exchange()
   only do this if rebalancing just occured
------------------------------------------------------------------------- */

void FixKineticsBalance::pre_neighbor()
{
  if (!pending) return;
  imbfinal = balance->imbalance_factor(maxloadperproc);
  pending = 0;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

inline double trim(double split, int n)
{
  int i = split * n;
  return (double)i / n;
}

void FixKineticsBalance::rebalance()
{
  imbprev = imbnow;

  // invoke balancer and reset comm->uniform flag

  int *sendproc;
  if (lbstyle == SHIFT) {
    itercount = balance->shift();
    comm->layout = LAYOUT_NONUNIFORM;
    for (int i = 0; i < comm->procgrid[0]; i++) {
      comm->xsplit[i] = trim(i * 1.0 / comm->procgrid[0], kinetics->nx);
    }
    for (int i = 0; i < comm->procgrid[1]; i++) {
      comm->ysplit[i] = trim(i * 1.0 / comm->procgrid[1], kinetics->ny);
    }
    for (int i = 0; i < comm->procgrid[2]; i++) {
      comm->zsplit[i] = trim(i * 1.0 / comm->procgrid[2], kinetics->nz);
    }
  } else if (lbstyle == BISECTION) {
    sendproc = balance->bisection();
    comm->layout = LAYOUT_TILED;
    int n[3];
    n[0] = kinetics->nx;
    n[1] = kinetics->ny;
    n[2] = kinetics->nz;
    double (*mysplit)[2] = comm->mysplit;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        mysplit[i][j] = trim(mysplit[i][j], n[i]);
      }
    }
    comm->rcbcutfrac = trim(comm->rcbcutfrac, n[comm->rcbcutdim]);
  }
  
  // output of new decomposition

  if (balance->outflag) balance->dumpout(update->ntimestep);

  // reset proc sub-domains
  // check and warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in exchange()

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();
  domain->subbox_too_small_check(neighbor->skin);

  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far
  // else allow caller's comm->exchange() to do it

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (wtflag) balance->fixstore->disable = 0;
  if (lbstyle == BISECTION) irregular->migrate_atoms(0,1,sendproc);
  else if (irregular->migrate_check()) irregular->migrate_atoms();
  if (wtflag) balance->fixstore->disable = 1;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // invoke KSpace setup_grid() to adjust to new proc sub-domains

  if (kspace_flag) force->kspace->setup_grid();

  // pending triggers pre_neighbor() to compute final imbalance factor
  // can only be done after atoms migrate in comm->exchange()

  pending = 1;

  kinetics->migrate();
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
------------------------------------------------------------------------- */

double FixKineticsBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
------------------------------------------------------------------------- */

double FixKineticsBalance::compute_vector(int i)
{
  if (i == 0) return maxloadperproc;
  if (i == 1) return (double) itercount;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixKineticsBalance::memory_usage()
{
  double bytes = irregular->memory_usage();
  if (balance->rcb) bytes += balance->rcb->memory_usage();
  return bytes;
}
