"""ASE LAMMPS Calculator Library Version"""

import os
import numpy as np
from numpy.linalg import norm
from lammps import lammps
from ase.calculators.calculator import Calculator
from ase.units import GPa
import ase.units
import ctypes

# TODO
# 1. should we make a new lammps object each time ?
# 2. upper triangular test does not look good
# 3. lmp object is not closed
# 4. need a routine to get the model back from lammps
# 5. if we send a command to lmps directly then the calculator does
#    not know about it and the energy could be wrong.

# 6. do we need a subroutine generator that converts a lammps string
#   into a python function that can be called


class LAMMPSlib(Calculator):
    r"""
    LAMMPSlib Interface Documentation
    
Introduction
============

LAMMPSlib is an interface and calculator for LAMMPS_. LAMMPSlib uses
the python interface that comes with LAMMPS to solve an atoms model
for energy, atom forces and cell stress. This calculator creates a
'.lmp' object which is a running lammps program, so further commands
can be sent to this object executed until it is explicitly closed. Any
additional variables calculated by lammps can also be extracted. This
is still experimental code.
    
Arguments
=========

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
                   appear in the atoms model

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


Requirements
============

To run this calculator you must have LAMMPS installed and compiled to
enable the python interface. See the LAMMPS manual.

If the following code runs then lammps is installed correctly.

   >>> from lammps import lammps
   >>> lmp = lammps()

The version of LAMMPS is also important. LAMMPSlib is suitable for
versions after approximately 2011. Prior to this the python interface
is slightly different from that used by LAMMPSlib. It is not difficult
to change to the earlier format.

LAMMPS and LAMMPSlib
====================

The LAMMPS calculator is another calculator that uses LAMMPS (the
program) to calculate the energy by generating input files and running
a separate LAMMPS job to perform the analysis. The output data is then
read back into python. LAMMPSlib makes direct use of the LAMMPS (the
program) python interface. As well as directly running any LAMMPS
comand line it allows the values of any of LAMMPS variables to be
extracted and returned to python.

Example
=======

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

    
Implementation
==============

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


Notes
=====

.. _LAMMPS: http://lammps.sandia.gov/

* Units: The default lammps_header sets the units to Angstrom and eV
  and for compatibility with ASE Stress is in GPa.

* The global energy is currently extracted from LAMMPS using
  extract_variable since lammps.lammps currently extract_global only
  accepts the following ['dt', 'boxxlo', 'boxxhi', 'boxylo', 'boxyhi',
  'boxzlo', 'boxzhi', 'natoms', 'nlocal'].

* If an error occurs while lammps is in control it will crash
  Python. Check the output of the log file to find the lammps error.

* If the are commands directly sent to the LAMMPS object this may
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
                       'atom_modify map array sort 0 0'])

    def set_cell(self, atoms, change=False):
        cell = self.convert_cell(atoms.get_cell())
        xhi = cell[0, 0]
        yhi = cell[1, 1]
        zhi = cell[2, 2]
        xy = cell[0, 1]
        xz = cell[0, 2]
        yz = cell[1, 2]

        if change:
            cell_cmd = 'change_box all     x final 0 {} y final 0 {} z final 0 {}      xy final {} xz final {} yz final {}'\
                .format(xhi, yhi, zhi, xy, xz, yz)
        else:
            # just in case we'll want to run with a funny shape box, and here command will only happen once, and before any calculation
            self.lmp.command('box tilt large')
            cell_cmd = 'region cell prism    0 {} 0 {} 0 {}     {} {} {}     units box'\
                .format(xhi, yhi, zhi, xy, xz, yz)

        self.lmp.command(cell_cmd)

    def calculate(self, atoms, properties, system_changes):
        self.propagate(atoms, properties, system_changes, 0)

    def propagate(self, atoms, properties, system_changes, n_steps, dt=None):

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

        self.atom_types = None
        self.coord_transform = None

	if not self.started:
            self.start_lammps()

	#NB
	if not self.initialized:
	   self.initialise_lammps(atoms)
        else: # still need to reset cell
           self.set_cell(atoms, change=True)

        pos = atoms.get_positions()

        # If necessary, transform the positions to new coordinate system
        if self.coord_transform != None:
            pos = np.dot(self.coord_transform , np.matrix.transpose(pos))
            pos = np.matrix.transpose(pos)

        # Convert ase position matrix to lammps-style position array
        lmp_positions = list(pos.ravel())

        # Convert that lammps-style array into a C object
        lmp_c_positions =\
            (ctypes.c_double * len(lmp_positions))(*lmp_positions)
#        self.lmp.put_coosrds(lmp_c_positions)
        self.lmp.scatter_atoms('x', 1, 3, lmp_c_positions)


        if n_steps > 0:
            vel = atoms.get_velocities()/(ase.units.Ang/(1.0e-12*ase.units.s))

            # If necessary, transform the velocities to new coordinate system
            if self.coord_transform != None:
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
            self.lmp.command('timestep %f' % (dt/(1.0e-12*ase.units.s)))
        self.lmp.command('run %d' % n_steps)

        if n_steps > 0:
            # TODO this must be slower than native copy, but why is it broken?
            pos = np.array([x for x in self.lmp.gather_atoms("x",1,3)]).reshape(-1,3)
            if self.coord_transform is not None:
                pos = np.dot(pos, self.coord_transform)
            atoms.set_positions(pos)
            vel = np.array([v for v in self.lmp.gather_atoms("v",1,3)]).reshape(-1,3)
            if self.coord_transform is not None:
                vel = np.dot(vel, self.coord_transform)
            atoms.set_velocities(vel*(ase.units.Ang/(1.0e-12*ase.units.s)))

        # Extract the forces and energy
#        if 'energy' in properties:
        self.results['energy'] = self.lmp.extract_variable('pe', None, 0)
#            self.results['energy'] = self.lmp.extract_global('pe', 0)
            
#        if 'stress' in properties:
        stress = np.empty(6)
        stress_vars = ['pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz']

        for i, var in enumerate(stress_vars):
            stress[i] = self.lmp.extract_variable(var, None, 0)
            
        # 1 bar (used by lammps for metal units) = 1e-4 GPa
        self.results['stress'] = stress * -1e-4 * GPa

#        if 'forces' in properties:
        f = np.zeros((len(atoms), 3))
        force_vars = ['fx', 'fy', 'fz']
        for i, var in enumerate(force_vars):
            f[:, i] = np.asarray(self.lmp.extract_variable(
                    var, 'all', 1)[:len(atoms)])
            
        self.results['forces'] = f

        if not self.parameters.keep_alive:
            self.lmp.close()

    def is_upper_triangular(self, mat):
        """test if 3x3 matrix is upper triangular"""
        
        def near0(x):
            """Test if a float is within .00001 of 0"""
            return abs(x) < 0.00001
        
        return near0(mat[1, 0]) and near0(mat[2, 0]) and near0(mat[2, 1])

    def convert_cell(self, ase_cell):
        """
        Convert a parallel piped (forming right hand basis)
        to lower triangular matrix LAMMPS can accept. This
        function transposes cell matrix so the bases are column vectors
        """
        cell = np.matrix.transpose(ase_cell)

        if not self.is_upper_triangular(cell):
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
            self.coord_transform = np.dot(tri_mat , trans)

            return tri_mat
        else:
            return cell

    def lammpsbc(self, pbc):
        if pbc:
            return 'p'
        else:
            return 's'

    def start_lammps(self):
        # start lammps process
        if self.parameters.log_file == None:
            cmd_args = ['-echo', 'log', '-log', 'none', '-screen', 'none']
        else:
            cmd_args = ['-echo', 'log', '-log', self.parameters.log_file,
                        '-screen', 'none']

        self.cmd_args = cmd_args

        if not hasattr(self, 'lmp'):
            self.lmp = lammps(self.parameters.lammps_name, self.cmd_args)

        # Use metal units: Angstrom, ps, and eV
        for cmd in self.parameters.lammps_header:
            self.lmp.command(cmd)

        self.started=True

    def initialise_lammps(self, atoms):

        # Initialising commands

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
        self.set_cell(atoms)

        # The default atom_types has atom type in alphabetic order
        # by atomic symbol
        symbols = np.asarray(atoms.get_chemical_symbols())

        # if the dictionary of types has not already been specified
        if self.atom_types == None:
            self.atom_types = {}
            atom_types = np.sort(np.unique(symbols))

            for i, sym in enumerate(atom_types):
                self.atom_types[sym] = i + 1

        # Initialize box
        n_types = len(self.atom_types)
        types_command = 'create_box {} cell'.format(n_types)
        self.lmp.command(types_command)

        # Initialize the atoms with their types
        # positions do not matter here
        self.lmp.command('echo none') # don't echo the atom positions
        for sym in symbols:
            cmd = 'create_atoms {} single 0.0 0.0 0.0  units box'.\
                format(self.atom_types[sym])
            self.lmp.command(cmd)

        self.lmp.command('echo log') # turn back on

        # Set masses
        masses = atoms.get_masses()
        for sym in self.atom_types:
            for i in range(len(atoms)):
                if symbols[i] == sym:
                    self.lmp.command('mass %d %f' % (self.atom_types[sym], masses[i]/(1.0e-3 * ase.units.kg /ase.units.mol))) # convert from amu (ASE) to g/mole (metal))
                    break

        # execute the user commands
        for cmd in self.parameters.lmpcmds:
            self.lmp.command(cmd)

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
        self.lmp.command('thermo_style custom pe pxx')
        
        self.lmp.command('variable fx atom fx')
        self.lmp.command('variable fy atom fy')
        self.lmp.command('variable fz atom fz')

        # do we need this if we extract from a global ?
        self.lmp.command('variable pe equal pe')

	self.initialized = True

#print('done loading lammpslib')
