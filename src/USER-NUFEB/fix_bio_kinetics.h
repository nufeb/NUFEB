/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics,FixKinetics)

#else

#ifndef SRC_FIX_KINETICS_H
#define SRC_FIX_KINETICS_H

#include "atom.h"
#include "bio.h"
#include "fix.h"
#include "decomp_grid.h"

namespace LAMMPS_NS {

class FixKinetics : public Fix, public DecompGrid<FixKinetics> {
  friend class DecompGrid<FixKinetics>;
  friend class FixKineticsEnergy;
  friend class FixKineticsMonod;
  friend class FixKineticsThermo;
  friend class FixKineticsDiffusion;
  friend class FixKineticsPH;
  friend class FixImmigration;
  friend class DumpBio;
  friend class FixKineticsBalance;

 public:
  FixKinetics(class LAMMPS *, int, char **);
  ~FixKinetics();
  int setmask();
  virtual void pre_force(int);
  void init();
  int modify_param(int, char **);
  void migrate();

  char **var;
  int *ivar;

  int nx, ny, nz, bnz;             // number of grids in x y z axis
  int bgrids;                      // # of non-boundary grids
  int ngrids;                      // # of grids

  double **nus;                    // nutrient concentration [nutrient][grid]
  double **nur;                    // nutrient consumption [nutrient][grid]
  double *nubs;                    // concentration in boundary layer [nutrient]
  double **fv;                     // velocity field [velo][grid]
  double **grid_yield;             // grid yield [type][grid]
  double ***activity;              // activities of chemical species [nutrient][5 charges][grid]
  double temp, rth;                // universal gas constant (thermodynamics) and temperature
  double **gibbs_cata;             // Gibbs free energy of catabolism [type][grid]
  double **gibbs_anab;             // Gibbs free energy of anabolism [type][grid]
  double **keq;                    // equilibrium constants [nutrient][4]
  double *sh;                      // concentration of hydrogen ion
  double **xdensity;               // grid biomass density [type][grid]; [0][grid] the overall density
  int *nuconv;                     // convergence flag
  double diff_dt;                  // diffusion timestep
  double blayer;                   // boundary layer
  double zhi,bzhi,zlo, xlo, xhi, ylo, yhi;
  double stepz, stepx, stepy;      // grid size
  int grow_flag;                   // microbe growth flag 1 = update biomass; 0 = solve reaction only, growth is negligible
  int demflag;                     // flag for DEM run
  double maxheight;                // maximum biofilm height
  int niter;                       // # of iterations
  int devery;                      // # of steps to call ph, thermo and form calculations

  int subn[3];                     // number of grids in x y axis for this proc
  int subnlo[3],subnhi[3];         // cell index of the subdomain lower and upper bound for each axis
  double sublo[3],subhi[3];        // subdomain lower and upper bound trimmed to the grid

  Grid<double, 3> grid;
  Subgrid<double, 3> subgrid;

  BIO *bio;
  class AtomVecBio *avec;
  class FixKineticsDiffusion *diffusion;
  class FixKineticsEnergy *energy;
  class FixKineticsMonod *monod;
  class FixKineticsPH *ph;
  class FixKineticsThermo *thermo;
  class FixFluid *nufebfoam;

  void init_param();
  void integration();
  void grow();
  double get_max_height();
  void update_bgrids();
  void update_xdensity();
  bool is_inside(int);
  int position(int);
  void reset_nur();
  void reset_isconv();

  Subgrid<double, 3> get_subgrid() const { return subgrid; }
  int get_elem_per_cell() const;
  template <typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnu; i++) {
	*result++ = nus[i][*it];
	*result++ = nur[i][*it];
      }
      if (energy) {
	for (int i = 1; i <= bio->nnu; i++) {
	  for (int j = 0; j < 5; j++) {
	    *result++ = activity[i][j][*it];
	  }
	}
	for (int i = 1; i <= atom->ntypes; i++) {
	  *result++ = grid_yield[i][*it];
	  *result++ = gibbs_cata[i][*it];
	  *result++ = gibbs_anab[i][*it];
	}
      }
      if (nufebfoam) {
	for (int i = 0; i < 3; i++) {
	  *result++ = fv[i][*it];
	}
      }
    }
    return result;
  }
  template <typename InputIterator0, typename InputIterator1>
  InputIterator1 unpack_cells(InputIterator0 first, InputIterator0 last, InputIterator1 input) {
    for (InputIterator0 it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnu; i++) {
	nus[i][*it] = *input++;
	nur[i][*it] = *input++;
      }
      if (energy) {
	for (int i = 1; i <= bio->nnu; i++) {
	  for (int j = 0; j < 5; j++) {
	    activity[i][j][*it] = *input++;
	  }
	}
	for (int i = 1; i <= atom->ntypes; i++) {
	  grid_yield[i][*it] = *input++;
	  gibbs_cata[i][*it] = *input++;
	  gibbs_anab[i][*it] = *input++;
	}
      }
      if (nufebfoam) {
	for (int i = 0; i < 3; i++) {
	  fv[i][*it] = *input++;
	}
      }
    }
    return input;
  }
  void resize(const Subgrid<double, 3> &);
};

}

#endif
#endif

