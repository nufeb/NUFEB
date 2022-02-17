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

#ifndef LMP_ERROR_H
#define LMP_ERROR_H

#include "pointers.h"

#ifdef LAMMPS_EXCEPTIONS
#include "exceptions.h"
#endif

namespace LAMMPS_NS {

class Error : protected Pointers {
 public:
  Error(class LAMMPS *);

  void universe_all(const char *, int, const char *);
  void universe_one(const char *, int, const char *);
  void universe_warn(const char *, int, const char *);

  void all(const char *, int, const char *);
  void one(const char *, int, const char *);
  void warning(const char *, int, const char *, int = 1);
  void message(const char *, int, const char *, int = 1);
  void done(int = 0); // 1 would be fully backwards compatible

#ifdef LAMMPS_EXCEPTIONS
  char *    get_last_error() const;
  ErrorType get_last_error_type() const;
  void   set_last_error(const char * msg, ErrorType type = ERROR_NORMAL);

 private:
  char * last_error_message;
  ErrorType last_error_type;
#endif
};

}

#endif

/* ERROR/WARNING messages:

*/
