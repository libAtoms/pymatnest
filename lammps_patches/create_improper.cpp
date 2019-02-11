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

#include <stdlib.h>
#include <string.h>
#include "create_improper.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "comm.h"
#include "group.h"
#include "special.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CreateImproper::CreateImproper(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateImproper::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_improper command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use create_improper unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use create_improper with non-molecular system");

  if (narg != 5) error->all(FLERR,"Illegal create_improper command");

  int *num_improper = atom->num_improper;
  int **improper_type = atom->improper_type;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  double newton_bond = force->newton_bond;
  int m;

  // parse args
  //

  int atype = force->inumeric(FLERR,arg[0]);
  tagint atom1 = force->tnumeric(FLERR,arg[1]);
  tagint atom2 = force->tnumeric(FLERR,arg[2]);
  tagint atom3 = force->tnumeric(FLERR,arg[3]);
  tagint atom4 = force->tnumeric(FLERR,arg[4]);


  if (atype <= 0 || atype > atom->nimpropertypes)
    error->all(FLERR,"Invalid improper type in create_improper command");

  // store state before improper creation

  bigint nimpropers_previous = atom->nimpropers;

  // check that atom ids are valid

if (atom1 <= 0 || atom1 > atom->map_tag_max ||
    atom2 <= 0 || atom2 > atom->map_tag_max ||
    atom3 <= 0 || atom3 > atom->map_tag_max )
  error->one(FLERR,"Invalid atom ID in create_improper command");
if (atype <= 0 || atype > atom->nimpropertypes)
  error->one(FLERR,"Invalid improper type in create_improper command");

if ((m = atom->map(atom2)) >= 0) {
    improper_type[m][num_improper[m]] = atype;
    improper_atom1[m][num_improper[m]] = atom1;
    improper_atom2[m][num_improper[m]] = atom2;
    improper_atom3[m][num_improper[m]] = atom3;
    improper_atom4[m][num_improper[m]] = atom4;
    num_improper[m]++;
}
if (newton_bond == 0) {
  if ((m = atom->map(atom1)) >= 0) {
      improper_type[m][num_improper[m]] = atype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
  }
  if ((m = atom->map(atom3)) >= 0) {
      improper_type[m][num_improper[m]] = atype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
  }

}

  // recount impropers

int nlocal = atom->nlocal;
  bigint nimpropers = 0;
  for (int i = 0; i < nlocal; i++) nimpropers += num_improper[i];

  MPI_Allreduce(&nimpropers,&atom->nimpropers,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (!force->newton_bond) atom->nimpropers /= 3;

  // print new improper count

  bigint nadd_impropers = atom->nimpropers - nimpropers_previous;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Added " BIGINT_FORMAT
              " impropers, new total = " BIGINT_FORMAT "\n",
              nadd_impropers,atom->nimpropers);
    }

    if (logfile) {
      fprintf(logfile,"Added " BIGINT_FORMAT
              " impropers, new total = " BIGINT_FORMAT "\n",
              nadd_impropers,atom->nimpropers);
    }
  }

  // re-trigger special list build

  Special special(lmp);
  special.build();
}
