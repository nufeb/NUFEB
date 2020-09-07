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

FixStyle(kinetics/diffusion,FixKineticsDiffusion)

#else

#ifndef SRC_FIX_KINETICS_DIFFUSIONS_H
#define SRC_FIX_KINETICS_DIFFUSIONS_H

#include "bio.h"
#include "fix.h"
#include "decomp_grid.h"

const int M2L = 1000;                   // convert from m3 to litre
const double MIN_NUS = 1e-20;           // minimal concentration

namespace LAMMPS_NS {
class AtomVecBio;
class BIO;
class FixKinetics;

class FixKineticsDiffusion: public Fix, public DecompGrid<FixKineticsDiffusion> {
  friend DecompGrid<FixKineticsDiffusion> ;

public:
  FixKineticsDiffusion(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion();

  char **var;
  int *ivar;

  double *rmass;
  int ntypes;
  int *mask;
  int *type;

  double stepx, stepy, stepz;             // grid size in x, y, z axis

  int xpbcflag, ypbcflag, zpbcflag;          // boundary condition flags for x, y,z. 0=PERIODIC 1=DIRICHLET
  int bulkflag;                           // 1=solve mass balance for bulk liquid
  int shearflag, dragflag, dcflag;        // flags for shear, drag(nufebfoam), and diffusion coefficent
  double dcratio;                         // ratio of diffusion coefficent in biomass region

  int *grid_type;                         // ghost grid flag [gird] 1=ghost gird, 0=non-ghost grid

  double srate;                           // shear rate
  double tol;                             // tolerance for convergence criteria

  double **nugrid;                        // nutrient concentration in ghost grid [nutrient][grid], unit in mol or kg/m3
  double **xgrid;                         // grid coordinate [gird][3]
  double **nuprev;                        // nutrient concentration in previous diffusion step
  double **nupenult;			  // nutrient concentration in n-2 step
  double **grid_diff_coeff;               // diffusion coeffs at each grid
  double vol;                             // grid volume

  double q, rvol, af;                     // parameters used for dynamic bulk
  int unit;                               // concentration unit 0=mol/l; 1=kg/m3

  int nx, ny, nz;                         // # of non-ghost grids in x, y and z axis
  int snxyz;                              // total # of local non-ghost grid
  int snxx, snyy, snzz;                   // # of local + ghost grids in x, y and z
  int snxx_yy;                            // # of local + ghost grids in xy
  int snxx_yy_zz;                         // total # of local + ghost grids

  double diff_dt;
  int close_flag;                        // flag for close system 1 = gradient is negligible (no diffusion)
  int* nuclose;
  int close_system;

  double xlo, xhi, ylo, yhi, zlo, zhi, bzhi;
  double init_nusbc;                     // inlet BC concentrations for dirichlet bc

  MPI_Request *requests;

  BIO *bio;
  FixKinetics *kinetics;
  AtomVecBio *avec;

  bool setup_exchange_flag; // flags that setup_exchange needs to be called in the next call to diffusion
  
  int setmask();
  void init();
  void init_setting();
  void init_grid();
  void closed_avg(int*, double);
  void closed_res(double, double);
  int *diffusion(int*, int, double);
  void update_nus();
  void update_grids();
  void update_diff_coeff();
  void compute_bc(double &, double *, int, double, int *);
  void compute_bulk();
  void compute_blayer();
  void compute_flux(double, double &, double *, double, int, int);
  void check_converge(int *);

  int get_index(int);
  void migrate(const Grid<double, 3> &, const Box<int, 3> &, const Box<int, 3> &);

  int get_elem_per_cell() const;
  template<typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnu; i++) {
        *result++ = nugrid[i][*it];
      }
    }
    return result;
  }
  template<typename InputIterator0, typename InputIterator1>
  InputIterator1 unpack_cells(InputIterator0 first, InputIterator0 last, InputIterator1 input) {
    for (InputIterator0 it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnu; i++) {
        nugrid[i][*it] = *input++;
      }
    }
    return input;
  }
  void resize(const Subgrid<double, 3> &);
};

}

#endif
#endif

