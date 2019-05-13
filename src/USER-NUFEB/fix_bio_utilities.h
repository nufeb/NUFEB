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

FixStyle(utilities,FixUtilities)

#else

#ifndef LMP_FIX_UTILITIES_H
#define LMP_FIX_UTILITIES_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixUtilities : public Fix {

 public:
  FixUtilities(class LAMMPS *, int, char **);
 ~FixUtilities();
  int setmask();
  void init();
  virtual void post_force(int);
  void get_floc (int, std::vector<int>&);
  void neighbor_list ();

 private:

  int nall;
  int *visit;
  double fourThirdsPI;
  std::vector<bigint> id;
  FILE* pFile;

  std::vector< std::vector<int> > list;

//  class NeighList *list;
};

}

#endif
#endif


