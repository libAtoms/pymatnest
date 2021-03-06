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

#ifdef COMMAND_CLASS

CommandStyle(create_improper,CreateImproper)

#else

#ifndef LMP_CREATE_IMPROPER_H
#define LMP_CREATE_IMPROPER_H

#include "pointers.h"

namespace LAMMPS_NS {

class CreateImproper : protected Pointers {
 public:
  CreateImproper(class LAMMPS *);
  void command(int, char **);

 private:
  inline int sbmask(int j) {
    return j >> SBBITS & 3;
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Create_improper command before simulation box is defined

Self-explanatory.

E: Cannot use create_improper unless atoms have IDs

This command requires a mapping from global atom IDs to local atoms,
but the atoms that have been defined have no IDs.

E: Cannot use create_improper with non-molecular system

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid improper type in create_improper command

Self-explanatory.

E: New improper exceeded impropers per atom in create_improper

See the read_data command for info on setting the "extra improper per
atom" header value to allow for additional impropers to be formed.

*/
