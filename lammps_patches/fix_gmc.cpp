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

#include <stdio.h>
#include <string.h>
#include "fix_gmc.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "compute.h"
#include "update.h"
#include "random_mars.h"
#include "memory.h"
#include "respa.h"
#include "error.h"
#include <math.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGMC::FixGMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"gmc") != 0 && narg < 5)
    error->all(FLERR,"Illegal fix gmc command");

  dynamic_group_allow = 1;
  time_integrate = 1;

  seed = utils::inumeric(FLERR,arg[3],true,lmp);
  Emax = utils::numeric(FLERR,arg[4],true,lmp);

  random = new RanMars(lmp,seed + comm->me);
  memory->create(dx,atom->nmax,3,"gmc:dx");


}

/* ---------------------------------------------------------------------- */

int FixGMC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGMC::init()
{

  int id = modify->find_compute("thermo_pe");

  modify->compute[id]->invoked_scalar = -1;
  pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep+1);


  int nlocal = atom->nlocal;
  double dx2sum = 0;

  for(int i = 0; i < nlocal; i++) {
    dx[i][0] = random->gaussian();
    dx[i][1] = random->gaussian();
    dx[i][2] = random->gaussian();
    dx2sum += dx[i][0]*dx[i][0]+dx[i][1]*dx[i][1]+dx[i][2]*dx[i][2];
  }

  dx2sum = sqrt(dx2sum);

  for(int i = 0; i < nlocal; i++) {
    dx[i][0] /= dx2sum;
    dx[i][1] /= dx2sum;
    dx[i][2] /= dx2sum;
  }

  step_size = update->dt;

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"fix gmc not compatible with RESPA");

}

void FixGMC::initial_integrate(int vflag)
{
  // update x of atoms in group

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
        

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      x[i][0] += step_size * dx[i][0];
      x[i][1] += step_size * dx[i][1];
      x[i][2] += step_size * dx[i][2];
    }

  int id = modify->find_compute("thermo_pe");

  modify->compute[id]->invoked_scalar = -1;
  pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep+1);


}

/* ---------------------------------------------------------------------- */

void FixGMC::final_integrate()
{

  // if potential energy is above Emax then want to modify dx with
  // forces to change trajectory
 
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double ecurrent = modify->compute[modify->find_compute("thermo_pe")]->compute_scalar();

  if ( ecurrent >= Emax) {
    double fsum = 0;
    double fhatdotdx = 0;
    for (int i = 0; i < nlocal; i++) 
        if (mask[i] & groupbit) {
          fsum += f[i][0]*f[i][0]+
                  f[i][1]*f[i][1]+
                  f[i][2]*f[i][2];
        }
    fsum = sqrt(fsum);
          
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
            fhatdotdx += f[i][0]/fsum*dx[i][0];
            fhatdotdx += f[i][1]/fsum*dx[i][1];
            fhatdotdx += f[i][2]/fsum*dx[i][2];
      }
    }

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
            dx[i][0] -= 2*f[i][0]/fsum*fhatdotdx;        
            dx[i][1] -= 2*f[i][1]/fsum*fhatdotdx;        
            dx[i][2] -= 2*f[i][2]/fsum*fhatdotdx;        
      }
    }
  }




}
FixGMC::~FixGMC()
{

  // delete temperature and pressure if fix created them
  delete random;

  memory->destroy(dx);

}


