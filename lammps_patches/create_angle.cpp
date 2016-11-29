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
#include "create_angle.h"
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

CreateAngle::CreateAngle(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateAngle::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_angle command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use create_angle unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use create_angle with non-molecular system");

  if (narg != 4) error->all(FLERR,"Illegal create_angle command");

  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  double newton_bond = force->newton_bond;
  int m;

  // parse args
  //

  int atype = force->inumeric(FLERR,arg[0]);
  tagint atom1 = force->tnumeric(FLERR,arg[1]);
  tagint atom2 = force->tnumeric(FLERR,arg[2]);
  tagint atom3 = force->tnumeric(FLERR,arg[3]);


  if (atype <= 0 || atype > atom->nangletypes)
    error->all(FLERR,"Invalid angle type in create_angle command");

  // store state before angle creation

  bigint nangles_previous = atom->nangles;

  // check that atom ids are valid

if (atom1 <= 0 || atom1 > atom->map_tag_max ||
    atom2 <= 0 || atom2 > atom->map_tag_max ||
    atom3 <= 0 || atom3 > atom->map_tag_max )
  error->one(FLERR,"Invalid atom ID in create_angle command");
if (atype <= 0 || atype > atom->nangletypes)
  error->one(FLERR,"Invalid angle type in create_angle command");

if ((m = atom->map(atom2)) >= 0) {
    angle_type[m][num_angle[m]] = atype;
    angle_atom1[m][num_angle[m]] = atom1;
    angle_atom2[m][num_angle[m]] = atom2;
    angle_atom3[m][num_angle[m]] = atom3;
    num_angle[m]++;
}
if (newton_bond == 0) {
  if ((m = atom->map(atom1)) >= 0) {
      angle_type[m][num_angle[m]] = atype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
  }
  if ((m = atom->map(atom3)) >= 0) {
      angle_type[m][num_angle[m]] = atype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
  }

}

  // recount angles

int nlocal = atom->nlocal;
  bigint nangles = 0;
  for (int i = 0; i < nlocal; i++) nangles += num_angle[i];

  MPI_Allreduce(&nangles,&atom->nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (!force->newton_bond) atom->nangles /= 3;

  // print new angle count

  bigint nadd_angles = atom->nangles - nangles_previous;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Added " BIGINT_FORMAT
              " angles, new total = " BIGINT_FORMAT "\n",
              nadd_angles,atom->nangles);
    }

    if (logfile) {
      fprintf(logfile,"Added " BIGINT_FORMAT
              " angles, new total = " BIGINT_FORMAT "\n",
              nadd_angles,atom->nangles);
    }
  }

  // re-trigger special list build

  Special special(lmp);
  special.build();
}
