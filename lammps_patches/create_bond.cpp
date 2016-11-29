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
#include "create_bond.h"
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

CreateBond::CreateBond(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateBond::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_bond command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use create_bond unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use create_bond with non-molecular system");

  if (narg != 3) error->all(FLERR,"Illegal create_bond command");

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  double newton_bond = force->newton_bond;
  int m;

  // parse args
  //

  int btype = force->inumeric(FLERR,arg[0]);
  tagint atom1 = force->tnumeric(FLERR,arg[1]);
  tagint atom2 = force->tnumeric(FLERR,arg[2]);

  if (btype <= 0 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in create_bond command");

  // store state before bond creation
  
  bigint nbonds_previous = atom->nbonds;

  // check that atom ids are valid

if (atom1 <= 0 || atom1 > atom->map_tag_max ||
    atom2 <= 0 || atom2 > atom->map_tag_max)
  error->one(FLERR,"Invalid atom ID in create_bond command");
if (btype <= 0 || btype > atom->nbondtypes)
  error->one(FLERR,"Invalid bond type in create_bond command");
//
if ((m = atom->map(atom1)) >= 0) {
    bond_type[m][num_bond[m]] = btype;
    bond_atom[m][num_bond[m]] = atom2;
    num_bond[m]++;
}
if (newton_bond == 0) {
  if ((m = atom->map(atom2)) >= 0) {
      bond_type[m][num_bond[m]] = btype;
      bond_atom[m][num_bond[m]] = atom1;
      num_bond[m]++;
  }
}

  // recount bonds

int nlocal = atom->nlocal;
  bigint nbonds = 0;
  for (int i = 0; i < nlocal; i++) nbonds += num_bond[i];

  MPI_Allreduce(&nbonds,&atom->nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (!force->newton_bond) atom->nbonds /= 2;

  // print new bond count

  bigint nadd_bonds = atom->nbonds - nbonds_previous;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Added " BIGINT_FORMAT
              " bonds, new total = " BIGINT_FORMAT "\n",
              nadd_bonds,atom->nbonds);
    }

    if (logfile) {
      fprintf(logfile,"Added " BIGINT_FORMAT
              " bonds, new total = " BIGINT_FORMAT "\n",
              nadd_bonds,atom->nbonds);
    }
  }

  // re-trigger special list build

  Special special(lmp);
  special.build();
}
