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
#include "create_dihedral.h"
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

CreateDihedral::CreateDihedral(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateDihedral::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_dihedral command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use create_dihedral unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use create_dihedral with non-molecular system");

  if (narg != 5) error->all(FLERR,"Illegal create_dihedral command");

  int *num_dihedral = atom->num_dihedral;
  int **dihedral_type = atom->dihedral_type;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  double newton_bond = force->newton_bond;
  int m;

  // parse args
  //

  int atype = force->inumeric(FLERR,arg[0]);
  tagint atom1 = force->tnumeric(FLERR,arg[1]);
  tagint atom2 = force->tnumeric(FLERR,arg[2]);
  tagint atom3 = force->tnumeric(FLERR,arg[3]);
  tagint atom4 = force->tnumeric(FLERR,arg[4]);


  if (atype <= 0 || atype > atom->ndihedraltypes)
    error->all(FLERR,"Invalid dihedral type in create_dihedral command");

  // store state before dihedral creation

  bigint ndihedrals_previous = atom->ndihedrals;

  // check that atom ids are valid

if (atom1 <= 0 || atom1 > atom->map_tag_max ||
    atom2 <= 0 || atom2 > atom->map_tag_max ||
    atom3 <= 0 || atom3 > atom->map_tag_max )
  error->one(FLERR,"Invalid atom ID in create_dihedral command");
if (atype <= 0 || atype > atom->ndihedraltypes)
  error->one(FLERR,"Invalid dihedral type in create_dihedral command");

if ((m = atom->map(atom2)) >= 0) {
    dihedral_type[m][num_dihedral[m]] = atype;
    dihedral_atom1[m][num_dihedral[m]] = atom1;
    dihedral_atom2[m][num_dihedral[m]] = atom2;
    dihedral_atom3[m][num_dihedral[m]] = atom3;
    dihedral_atom4[m][num_dihedral[m]] = atom4;
    num_dihedral[m]++;
}
if (newton_bond == 0) {
  if ((m = atom->map(atom1)) >= 0) {
      dihedral_type[m][num_dihedral[m]] = atype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
  }
  if ((m = atom->map(atom3)) >= 0) {
      dihedral_type[m][num_dihedral[m]] = atype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
  }

}

  // recount dihedrals

int nlocal = atom->nlocal;
  bigint ndihedrals = 0;
  for (int i = 0; i < nlocal; i++) ndihedrals += num_dihedral[i];

  MPI_Allreduce(&ndihedrals,&atom->ndihedrals,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (!force->newton_bond) atom->ndihedrals /= 3;

  // print new dihedral count

  bigint nadd_dihedrals = atom->ndihedrals - ndihedrals_previous;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Added " BIGINT_FORMAT
              " dihedrals, new total = " BIGINT_FORMAT "\n",
              nadd_dihedrals,atom->ndihedrals);
    }

    if (logfile) {
      fprintf(logfile,"Added " BIGINT_FORMAT
              " dihedrals, new total = " BIGINT_FORMAT "\n",
              nadd_dihedrals,atom->ndihedrals);
    }
  }

  // re-trigger special list build

  Special special(lmp);
  special.build();
}
