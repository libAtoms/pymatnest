"""ASE LAMMPS Calculator Library Version"""

import os
import ctypes
import operator

import numpy as np
from numpy.linalg import norm

from lammps import lammps
from ase.calculators.calculator import Calculator
from ase.data import chemical_symbols, atomic_masses
import ase.units
import re

# TODO
# 1. should we make a new lammps object each time ?
# 2. upper triangular test does not look good
# 3. lmp object is not closed
# 4. need a routine to get the model back from lammps
# 5. if we send a command to lmps directly then the calculator does
#    not know about it and the energy could be wrong.

# 6. do we need a subroutine generator that converts a lammps string
#   into a python function that can be called

def is_upper_triangular(mat):
    """test if 3x3 matrix is upper triangular"""

    def near0(x):
        """Test if a float is within .00001 of 0"""
        return abs(x) < 0.00001

    return near0(mat[1, 0]) and near0(mat[2, 0]) and near0(mat[2, 1])

def convert_cell(ase_cell):
    """
    Convert a parallel piped (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """
    cell = np.matrix.transpose(ase_cell)

    if not is_upper_triangular(cell):
        # rotate bases into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = cell[:, 0]
        B = cell[:, 1]
        C = cell[:, 2]
        tri_mat[0, 0] = norm(A)
        Ahat = A / norm(A)
        AxBhat = np.cross(A, B) / norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(ase_cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans = trans / volume
        coord_transform = np.dot(tri_mat , trans)

        return tri_mat, coord_transform
    else:
        return cell, None


lammps_real = {
      "mass" : 0.001 * ase.units.kg / ase.units.mol,
      "distance" : ase.units.Angstrom,
      "time" : ase.units.fs,
      "energy" : ase.units.kcal/ase.units.mol,
      "velocity": ase.units.Angstrom / ase.units.fs,
      "force": ase.units.kcal/ase.units.mol/ase.units.Angstrom,
      "pressure" : 101325 * ase.units.Pascal
      }

lammps_metal = {
      "mass" : 0.001 * ase.units.kg / ase.units.mol,
      "distance" : ase.units.Angstrom,
      "time" : 1e-12 * ase.units.second,
      "energy" : ase.units.eV,
      "velocity": ase.units.Angstrom / (1e-12*ase.units.second),
      "force": ase.units.eV/ase.units.Angstrom,
      "pressure" : 1e5 * ase.units.Pascal
      }

lammps_units={"real":lammps_real,
              "metal":lammps_metal}

def unit_convert(quantity, units='metal'):
   try:
      return lammps_units[units][quantity]
   except:
      raise NotImplementedError("Unit {} in unit system {} is not implemented.".format(quantity,units))

class LAMMPSlib(Calculator):

    r"""
    LAMMPSlib Interface Documentation

**Introduction**

LAMMPSlib is an interface and calculator for LAMMPS_. LAMMPSlib uses
the python interface that comes with LAMMPS to solve an atoms model
for energy, atom forces and cell stress. This calculator creates a
'.lmp' object which is a running lammps program, so further commands
can be sent to this object executed until it is explicitly closed. Any
additional variables calculated by lammps can also be extracted. This
is still experimental code.

**Arguments**

=================  ==========================================================
Keyword                               Description
=================  ==========================================================
``lmpcmds``        list of strings of LAMMPS commands. You need to supply
                   enough to define the potential to be used e.g.

                   ["pair_style eam/alloy",
                    "pair_coeff * * potentials/NiAlH_jea.eam.alloy Ni Al"]

``atom_types``     dictionary of "atomic_symbol":lammps_atom_type pairs,
                   e.g. {'Cu':1} to bind copper to lammps atom type 1.
                   Default method assigns lammps atom types in order that they
                   appear in the atoms model. Mandatory.

``log_file``       string
                   path to the desired LAMMPS log file

``lammps_header``  string to use for lammps setup. Default is to use
                   metal units and simple atom simulation.

                   lammps_header=['units metal',
                                  'atom_style atomic',
                                  'atom_modify map array sort 0 0'])

``keep_alive``     Boolean
                   whether to keep the lammps routine alive for more commands

=================  ==========================================================


**Requirements**

To run this calculator you must have LAMMPS installed and compiled to
enable the python interface. See the LAMMPS manual.

If the following code runs then lammps is installed correctly.

   >>> from lammps import lammps
   >>> lmp = lammps()

The version of LAMMPS is also important. LAMMPSlib is suitable for
versions after approximately 2011. Prior to this the python interface
is slightly different from that used by LAMMPSlib. It is not difficult
to change to the earlier format.

**LAMMPS and LAMMPSlib**

The LAMMPS calculator is another calculator that uses LAMMPS (the
program) to calculate the energy by generating input files and running
a separate LAMMPS job to perform the analysis. The output data is then
read back into python. LAMMPSlib makes direct use of the LAMMPS (the
program) python interface. As well as directly running any LAMMPS
comand line it allows the values of any of LAMMPS variables to be
extracted and returned to python.

**Example**

::

    from ase import Atom, Atoms
    from lammpslib import LAMMPSlib

    cmds = ["pair_style eam/alloy",
            "pair_coeff * * NiAlH_jea.eam.alloy Al H"]

    a = 4.05
    al = Atoms([Atom('Al')], cell=(a, a, a), pbc=True)
    h = Atom([Atom('H')])
    alh = al + h

    lammps = LAMMPSlib(lmpcmds = cmds, logfile='test.log')

    alh.set_calculator(lammps)
    print "Energy ", alh.get_potential_energy()


**Implementation**

LAMMPS provides a set of python functions to allow execution of the
underlying C++ LAMMPS code. The functions used by the LAMMPSlib
interface are::

    from lammps import lammps

    lmp = lammps(cmd_args) # initiate LAMMPS object with command line args

    lmp.scatter_atoms('x',1,3,positions) # atom coords to LAMMPS C array
    lmp.command(cmd) # executes a one line cmd string
    lmp.extract_variable(...) # extracts a per atom variable
    lmp.extract_global(...) # extracts a global variable
    lmp.close() # close the lammps object

For a single atom model the following lammps file commands would be run
by invoking the get_potential_energy() method::

    units metal
    atom_style atomic
    atom_modify map array sort 0 0

    region cell prism 0 xhi 0 yhi 0 zhi xy xz yz units box
    create_box 1 cell
    create_atoms 1 single 0 0 0 units box
    mass * 1.0

    ## user lmpcmds get executed here
    pair_style eam/alloy
    pair_coeff * * lammps/potentials/NiAlH_jea.eam.alloy Al
    ## end of user lmmpcmds

    run 0


**Notes**

.. _LAMMPS: http://lammps.sandia.gov/

* Units: The default lammps_header sets the units to Angstrom and eV
  and for compatibility with ASE Stress is in GPa.

* The global energy is currently extracted from LAMMPS using
  extract_variable since lammps.lammps currently extract_global only
  accepts the following ['dt', 'boxxlo', 'boxxhi', 'boxylo', 'boxyhi',
  'boxzlo', 'boxzhi', 'natoms', 'nlocal'].

* If an error occurs while lammps is in control it will crash
  Python. Check the output of the log file to find the lammps error.

* If the are commands direfctly sent to the LAMMPS object this may
  change the energy value of the model. However the calculator will not
  know of it and still return the original energy value.

End LAMMPSlib Interface Documentation

    """

    implemented_properties = ['energy', 'forces', 'stress']

    #NB
    started = False
    initialized = False

    default_parameters = dict(
        atom_types=None,
        log_file=None,
        lammps_name='',
        keep_alive=False,
        lammps_header=['units metal',
                       'atom_style atomic',
                       'atom_modify map array sort 0 0'],
        boundary=True,
        create_box=True,
        create_atoms=True,
        read_molecular_info=False,
        comm=None)

    def set_bonds(self, atoms):

        lbox = np.amin([atoms.get_cell()[0,0], atoms.get_cell()[1,1], atoms.get_cell()[2,2]])
        max_bond_length = lbox/2
        self.lmp.command('neighbor {} bin'.format(max_bond_length));

        pos = atoms.get_positions() / unit_convert("distance", self.units)
        for i in range(len(atoms)):
            if atoms.arrays['bonds'][i] != '_':
                for j in atoms.arrays['bonds'][i].split(','):
                    m = re.match('(\d+)\((\d+)\)',j)
                    distance = norm(pos[i,:]-pos[int(m.group(1)),:])
                    self.lmp.command('group g1 id {} '.format(i+1))
                    self.lmp.command('group g2 id {} '.format(int(m.group(1))+1))
                    self.lmp.command('create_bonds g1 g2 {} 0.0 {} '.format(int(m.group(2)),distance+1))
                    self.lmp.command('group g1 delete')
                    self.lmp.command('group g2 delete')

    def set_cell(self, atoms, change=False):
        lammps_cell, self.coord_transform = convert_cell(atoms.get_cell())
        xhi = lammps_cell[0, 0]
        yhi = lammps_cell[1, 1]
        zhi = lammps_cell[2, 2]
        xy = lammps_cell[0, 1]
        xz = lammps_cell[0, 2]
        yz = lammps_cell[1, 2]

        if change:
            cell_cmd = 'change_box all     x final 0 {} y final 0 {} z final 0 {}      xy final {} xz final {} yz final {}'\
                .format(xhi, yhi, zhi, xy, xz, yz)
        else:
            # just in case we'll want to run with a funny shape box, and here command will only happen once, and before any calculation
            if self.parameters.create_box:
                self.lmp.command('box tilt large')
            cell_cmd = 'region cell prism    0 {} 0 {} 0 {}     {} {} {}     units box'\
                .format(xhi, yhi, zhi, xy, xz, yz)

        self.lmp.command(cell_cmd)

    def set_lammps_pos(self, atoms):
        pos = atoms.get_positions() / unit_convert("distance", self.units)

        # If necessary, transform the positions to new coordinate system
        if self.coord_transform is not None:
            pos = np.dot(self.coord_transform , np.matrix.transpose(pos))
            pos = np.matrix.transpose(pos)

        # Convert ase position matrix to lammps-style position array
        lmp_positions = list(pos.ravel())

        # Convert that lammps-style array into a C object
        lmp_c_positions =\
            (ctypes.c_double * len(lmp_positions))(*lmp_positions)
#        self.lmp.put_coosrds(lmp_c_positions)
        self.lmp.scatter_atoms('x', 1, 3, lmp_c_positions)


    def calculate(self, atoms, properties, system_changes):
        self.propagate(atoms, properties, system_changes, 0)

    def propagate(self, atoms, properties, system_changes, n_steps, dt=None, dt_not_real_time=False):

        """"atoms: Atoms object
            Contains positions, unit-cell, ...
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these five: 'positions', 'numbers', 'cell',
            'pbc', 'charges' and 'magmoms'.
        """
        if len(system_changes) == 0:
            return

        self.coord_transform = None

        if not self.started:
            self.start_lammps()

        ####################################################################################################
        #NB
        if not self.initialized:
            self.initialise_lammps(atoms)
        else: # still need to reset cell
            # reset positions so that if they are cray from last propagation, change_box (in set_cell()) won't hang
            # could do this only after testing for crazy positions?
            # could also use scatter_atoms() to set values (requires MPI comm), or extra_atoms() to get pointers to local 
            #     data structures to zero, but then will have to be careful with parallelism
            self.lmp.command("set atom * x 0.0 y 0.0 z 0.0")
            self.set_cell(atoms, change=True)

        if self.parameters.atom_types is None:
           raise NameError("atom_types are mandatory.")

        do_rebuild = False
        do_redo_atom_types = False
        try:
            do_rebuild = ( len(atoms.numbers) != len(self.previous_atoms_numbers) ) or ( "numbers" in system_changes )
            if not do_rebuild:
                do_redo_atom_types = ( atoms.numbers != self.previous_atoms_numbers ).any()
        except Exception:
           pass

        self.lmp.command('echo none') # don't echo the atom positions
        if do_rebuild:
           self.rebuild(atoms)
        elif do_redo_atom_types:
           self.redo_atom_types(atoms)
        self.lmp.command('echo log') # switch back log

        self.set_lammps_pos(atoms)

        if n_steps > 0:
            vel = atoms.get_velocities() / unit_convert("velocity", self.units)

            # If necessary, transform the velocities to new coordinate system
            if self.coord_transform is not None:
                vel = np.dot(self.coord_transform , np.matrix.transpose(vel) )
                vel = np.matrix.transpose(vel)

            # Convert ase velocities matrix to lammps-style velocities array
            lmp_velocities = list(vel.ravel())

            # Convert that lammps-style array into a C object
            lmp_c_velocities =\
                (ctypes.c_double * len(lmp_velocities))(*lmp_velocities)
#            self.lmp.put_coosrds(lmp_c_velocities)
            self.lmp.scatter_atoms('v', 1, 3, lmp_c_velocities)

        # Run for 0 time to calculate
        if dt is not None:
            if dt_not_real_time:
                self.lmp.command('timestep %.30f' % dt)
            else:
                self.lmp.command('timestep %.30f' % ( dt/unit_convert("time", self.units)) )
        self.lmp.command('run %d' % n_steps)

        if n_steps > 0:
            # TODO this must be slower than native copy, but why is it broken?
            pos = np.array([x for x in self.lmp.gather_atoms("x",1,3)]).reshape(-1,3)
            if self.coord_transform is not None:
                pos = np.dot(pos, self.coord_transform)
            atoms.set_positions(pos * unit_convert("distance", self.units))
            vel = np.array([v for v in self.lmp.gather_atoms("v",1,3)]).reshape(-1,3)
            if self.coord_transform is not None:
                vel = np.dot(vel, self.coord_transform)
            atoms.set_velocities(vel * unit_convert("velocity", self.units))

        # Extract the forces and energy
#        if 'energy' in properties:
        self.results['energy'] = self.lmp.extract_variable('pe', None, 0) * unit_convert("energy", self.units)
#            self.results['energy'] = self.lmp.extract_global('pe', 0)

#        if 'stress' in properties:
        stress = np.empty(6)
        # stress_vars = ['pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz']
        stress_vars = ['pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy']

        for i, var in enumerate(stress_vars):
            stress[i] = self.lmp.extract_variable(var, None, 0)

        stress_mat = np.zeros( (3,3) )
        stress_mat[0,0] = stress[0]
        stress_mat[1,1] = stress[1]
        stress_mat[2,2] = stress[2]
        stress_mat[1,2] = stress[3]
        stress_mat[2,1] = stress[3]
        stress_mat[0,2] = stress[4]
        stress_mat[2,0] = stress[4]
        stress_mat[0,1] = stress[5]
        stress_mat[1,0] = stress[5]
        if self.coord_transform is not None:
            stress_mat = np.dot(self.coord_transform.T, np.dot(stress_mat, self.coord_transform))
        stress[0] = stress_mat[0,0]
        stress[1] = stress_mat[1,1]
        stress[2] = stress_mat[2,2]
        stress[3] = stress_mat[1,2]
        stress[4] = stress_mat[0,2]
        stress[5] = stress_mat[0,1]

        self.results['stress'] = stress * (-unit_convert("pressure", self.units))

#        if 'forces' in properties:
        f = np.zeros((len(atoms), 3))
        force_vars = ['fx', 'fy', 'fz']
        for i, var in enumerate(force_vars):
            f[:, i] = np.asarray(self.lmp.extract_variable(
                    var, 'all', 1)[:len(atoms)]) * unit_convert("force", self.units)

        if self.coord_transform is not None:
            self.results['forces'] = np.dot(f, self.coord_transform)
        else:
            self.results['forces'] = f.copy()

        if not self.parameters.keep_alive:
            self.lmp.close()

    def lammpsbc(self, pbc):
        if pbc:
            return 'p'
        else:
            return 's'

    def rebuild(self,atoms):

       try:
          n_diff = len(atoms.numbers) - len(self.previous_atoms_numbers)
       except:
          n_diff = len(atoms.numbers)

       if n_diff > 0:
          if any([("reax/c" in cmd) for cmd in self.parameters.lmpcmds]):
             self.lmp.command("pair_style lj/cut 2.5")
             self.lmp.command("pair_coeff * * 1 1")

             for cmd in self.parameters.lmpcmds:
                if ("pair_style" in cmd) or ("pair_coeff" in cmd):
                   self.lmp.command(cmd)

          cmd = "create_atoms 1 random {} 1 NULL".format(n_diff)
          self.lmp.command(cmd)
       elif n_diff < 0:
          cmd = "group delatoms id {}:{}".format(len(atoms.numbers)+1,len(self.previous_atoms_numbers))
          self.lmp.command(cmd)
          cmd = "delete_atoms group delatoms"
          self.lmp.command(cmd)

       self.redo_atom_types(atoms)

    def redo_atom_types(self,atoms):

       current_types = { (i+1,self.parameters.atom_types[sym]) for i,sym in enumerate( atoms.get_chemical_symbols() ) }

       try:
          previous_types = { (i+1,self.parameters.atom_types[ chemical_symbols[Z] ])
                              for i,Z in enumerate( self.previous_atoms_numbers ) }
       except:
          previous_types = set()

       for (i,i_type) in current_types - previous_types:
          cmd = "set atom {} type {}".format(i,i_type)
          self.lmp.command(cmd)

       self.previous_atoms_numbers = atoms.numbers.copy()

    def restart_lammps(self, atoms):
        if self.started:
            self.lmp.command("clear")
        # hope there's no other state to be reset
        self.started=False
        self.initialized=False
        self.previous_atoms_numbers = []
        self.start_lammps()
        self.initialise_lammps(atoms)

    def start_lammps(self):
        # start lammps process
        if self.parameters.log_file is None:
            cmd_args = ['-echo', 'log', '-log', 'none', '-screen', 'none', '-nocite']
        else:
            cmd_args = ['-echo', 'log', '-log', self.parameters.log_file,
                        '-screen', 'none','-nocite']

        self.cmd_args = cmd_args

        if not hasattr(self, 'lmp'):
            self.lmp = lammps(self.parameters.lammps_name, self.cmd_args, comm=self.parameters.comm)

        # Use metal units: Angstrom, ps, and eV
        for cmd in self.parameters.lammps_header:
            self.lmp.command(cmd)

        for cmd in self.parameters.lammps_header:
           if "units" in cmd:
              self.units = cmd.split()[1]

        if hasattr(self.parameters, "lammps_header_extra") and self.parameters.lammps_header_extra is not None:
            for cmd in self.parameters.lammps_header_extra:
                self.lmp.command(cmd)

        self.started=True

    def initialise_lammps(self, atoms):

        # Initialising commands
        if self.parameters.boundary:
            # if the boundary command is in the supplied commands use that
            # otherwise use atoms pbc
            pbc = atoms.get_pbc()
            for cmd in self.parameters.lmpcmds:
                if 'boundary' in cmd:
                    break
            else:
                self.lmp.command(
                    'boundary ' + ' '.join([self.lammpsbc(bc) for bc in pbc]))

        # Initialize cell
        self.set_cell(atoms, change=not self.parameters.create_box)

        if self.parameters.atom_types  is None:
           raise NameError("atom_types are mandatory.")

        # Collect chemical symbols
        symbols = np.asarray(atoms.get_chemical_symbols())

        # Initialize box
        if self.parameters.create_box:
           n_types = len(self.parameters.atom_types)
           types_command = 'create_box {} cell'.format(n_types)
           # Hard-coded some extra memory for different bond and angle types, and bonds and angles for each particle
           if self.parameters.read_molecular_info and 'angles' in atoms.arrays:
               types_command += ' angle/types 1 extra/angle/per/atom 2'
           if self.parameters.read_molecular_info and 'bonds' in atoms.arrays:
               types_command += ' bond/types 1 extra/bond/per/atom 4'
           self.lmp.command(types_command)

        # Initialize the atoms with their types
        # positions do not matter here
        if self.parameters.create_atoms:
           self.lmp.command('echo none') # don't echo the atom positions
           self.rebuild(atoms)
           self.lmp.command('echo log') # turn back on

        # execute the user commands
        for cmd in self.parameters.lmpcmds:
            self.lmp.command(cmd)

        # Set masses after user commands, to override EAM provided masses, e.g.
        masses = atoms.get_masses()
        for sym in self.parameters.atom_types:
            for i in range(len(atoms)):
                if symbols[i] == sym:
                    # convert from amu (ASE) to lammps mass unit)
                    self.lmp.command('mass %d %.30f' % (self.parameters.atom_types[sym], masses[i] /
                                                     unit_convert("mass", self.units) ))
                    break

        # Define force & energy variables for extraction
        self.lmp.command('variable pxx equal pxx')
        self.lmp.command('variable pyy equal pyy')
        self.lmp.command('variable pzz equal pzz')
        self.lmp.command('variable pxy equal pxy')
        self.lmp.command('variable pxz equal pxz')
        self.lmp.command('variable pyz equal pyz')

        # I am not sure why we need this next line but LAMMPS will
        # raise an error if it is not there. Perhaps it is needed to
        # ensure the cell stresses are calculated
        self.lmp.command('thermo_style custom pe pxx emol')
        self.lmp.command('thermo 1 ')

        self.lmp.command('variable fx atom fx')
        self.lmp.command('variable fy atom fy')
        self.lmp.command('variable fz atom fz')

        # do we need this if we extract from a global ?
        self.lmp.command('variable pe equal pe')

        self.lmp.command("neigh_modify delay 0 every 1 check yes")

        # read in bonds if there are bonds from the ase-atoms object if the atoms flag is set
        if self.parameters.read_molecular_info and 'bonds' in atoms.arrays:
            # Found we needed to have atom positions set for LAMMPS to properly bond atoms
            self.set_lammps_pos(atoms)
            self.set_bonds(atoms)

        self.initialized = True

#print('done loading lammpslib')

def write_lammps_data(filename, atoms, atom_types, comment=None, cutoff=None,
                      molecule_ids=None, charges=None, units='metal',
                      bond_types=None, angle_types=None, dihedral_types=None):

    if isinstance(filename, basestring):
        fh = open(filename, 'w')
    else:
        fh = filename

    if comment is None:
        comment = 'lammpslib autogenerated data file'
    fh.write(comment.strip() + '\n\n')

    fh.write('{0} atoms\n'.format(len(atoms)))
    fh.write('{0} atom types\n'.format(len(atom_types)))

    if bond_types:
        from matscipy.neighbours import neighbour_list
        i_list, j_list = neighbour_list('ij', atoms, cutoff)
        print 'Bonds:'
        bonds = []
        for bond_type, (Z1, Z2) in enumerate(bond_types):
            bond_mask = (atoms.numbers[i_list] == Z1) & (atoms.numbers[j_list] == Z2)
            print (Z1, Z2), bond_mask.sum()
            for (I, J) in zip(i_list[bond_mask], j_list[bond_mask]):
                #NB: LAMMPS uses 1-based indices for bond types and particle indices
                bond = (bond_type+1, I+1, J+1)
                bonds.append(bond)
        print
        if len(bonds) > 0:
            fh.write('{0} bonds\n'.format(len(bonds)))
            fh.write('{0} bond types\n'.format(len(bond_types)))

    if angle_types:
        print 'Angles:'
        angle_count = { angle : 0 for angle in angle_types }
        angles = []
        for I in range(len(atoms)):
            for J in j_list[i_list == I]:
                for K in j_list[i_list == J]:
                    if J < K:
                        continue
                    Zi, Zj, Zk = atoms.numbers[[I, J, K]]
                    if (Zj, Zi, Zk) in angle_types:
                        angle = (angle_types.index((Zj, Zi, Zk))+1, J+1, I+1, K+1)
                        angle_count[(Zj, Zi, Zk)] += 1
                        angles.append(angle)
        for angle in angle_types:
            print angle, angle_count[angle]
        print
        if len(angles) > 0:
            fh.write('{0} angles\n'.format(len(angles)))
            fh.write('{0} angle types\n'.format(len(angle_types)))

    if dihedral_types:
        print 'Dihedrals:'
        dihedral_count = { dihedral : 0 for dihedral in dihedral_types }
        dihedrals = []
        for I in range(len(atoms)):
            for J in j_list[i_list == I]:
                for K in j_list[i_list == J]:
                    for L in j_list[i_list == K]:
                        Zi, Zj, Zk, Zl = atoms.numbers[[I, J, K, L]]
                        if (Zi, Zj, Zk, Zl) in dihedral_types:
                            dihedral = (dihedral_types.index((Zi, Zj, Zk, Zl))+1,
                                        I+1, J+1, K+1, L+1)
                            dihedral_count[(Zi, Zj, Zk, Zl)] += 1
                            dihedrals.append(dihedral)
        for dihedral in dihedral_types:
            print dihedral, dihedral_count[dihedral]
        print
        if len(dihedrals) > 0:
            fh.write('{0} dihedrals\n'.format(len(dihedrals)))
            fh.write('{0} dihedral types\n'.format(len(dihedral_types)))

    fh.write('\n')
    cell, coord_transform = convert_cell(atoms.get_cell())
    fh.write('{0:16.8e} {1:16.8e} xlo xhi\n'.format(0.0, cell[0, 0]))
    fh.write('{0:16.8e} {1:16.8e} ylo yhi\n'.format(0.0, cell[1, 1]))
    fh.write('{0:16.8e} {1:16.8e} zlo zhi\n'.format(0.0, cell[2, 2]))
    fh.write('{0:16.8e} {1:16.8e} {2:16.8e} xy xz yz\n'.format(cell[0, 1], cell[0, 2], cell[1, 2]))

    fh.write('\nMasses\n\n')
    sym_mass = {}
    masses = atoms.get_masses()
    symbols = atoms.get_chemical_symbols()
    for sym in atom_types:
        for i in range(len(atoms)):
            if symbols[i] == sym:
                sym_mass[sym] = masses[i] / unit_convert("mass", units)
                break
            else:
                sym_mass[sym] = atomic_masses[chemical_symbols.index(sym)] / unit_convert("mass", units)

    for (sym, typ) in sorted(atom_types.items(), key=operator.itemgetter(1)):
        fh.write('{0} {1}\n'.format(typ, sym_mass[sym]))

    fh.write('\nAtoms # full\n\n')
    if molecule_ids is None:
        molecule_ids = np.zeros(len(atoms), dtype=int)
    if charges is None:
        charges = atoms.get_initial_charges()
    for i, (sym, mol, q, pos) in enumerate(zip(symbols, molecule_ids,
                                               charges, atoms.get_positions())):
        typ = atom_types[sym]
        fh.write('{0} {1} {2} {3:16.8e} {4:16.8e} {5:16.8e} {6:16.8e}\n'
                 .format(i+1, mol, typ, q, pos[0], pos[1], pos[2]))

    if bond_types and len(bonds) > 0:
        fh.write('\nBonds\n\n')
        for idx, bond in enumerate(bonds):
            fh.write('{0} {1} {2} {3}\n'
                     .format(*[idx+1] + list(bond)))

    if angle_types and len(angles) > 0:
        fh.write('\nAngles\n\n')
        for idx, angle in enumerate(angles):
            fh.write('{0} {1} {2} {3} {4}\n'
                     .format(*[idx+1] + list(angle)))

    if dihedral_types and len(dihedrals) > 0:
        fh.write('\nDihedrals\n\n')
        for idx, dihedral in enumerate(dihedrals):
            fh.write('{0} {1} {2} {3} {4} {5}\n'
                     .format(*[idx+1] + list(dihedral)))

    if isinstance(filename, basestring):
        fh.close()


