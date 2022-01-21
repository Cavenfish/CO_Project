import sys, time
import numpy as np
from ase.io import read, write
from ase import units, Atoms
from ase.optimize import BFGS
from .MvH_CO_JM8 import MvH_CO
from ase.visualize import view
from ase.vibrations import Vibrations
from ase.md.verlet import VelocityVerlet
from ase.io.trajectory import Trajectory

def view_xyz(xyz):
    system = read(xyz)
    view(system)
    return 

def prep_system(xyz):
    system = read(xyz)
    calc   = MvH_CO(atoms=system)
    system.set_calculator(calc)
    return system, calc

def geo_opt(xyz):
    #Read in system and set van Hemert calculator
    system, calc = prep_system(xyz)

    #Make trajectory string
    traj  = xyz.replace('.xyz', '_opt.traj')
    
    #Run BFGS optimization of geometry
    opt   = BFGS(system, trajectory=traj)
    opt.run(fmax=0.0001)

    #Make XYZ file of optimized system
    traj  = Trajectory(traj)
    atoms = traj[-1]
    opt_f = xyz.replace('.xyz', '_opt.xyz')
    write(opt_f, atoms)

    return

def calc_vibs(xyz):
    #Read in system and set van Hemert calculator
    system, calc = prep_system(xyz)

    #Run vibrational analysis
    vib = Vibrations(system, delta=0.0001)
    vib.run()
    vib.summary()

    return

def get_system_properties(xyz):
    #Read in system and set van Hemert calculator
    system, calc = prep_system(xyz)

    #Get system potential energy and track time of calculation
    b4          = time.time()
    E_pot       = system.get_potential_energy()
    E_pot_time  = time.time() - b4

    #Pretty Header
    print('\n\n---Starting Printout---\n\n')

    #Pretty print potential energy
    print(' Potential Energy: %.4f\n' % E_pot,
           'Runtime         : %.1f (s)\n\n' % E_pot_time)

    #Get energy contributions
    E_intra, E_pair, E_ex, E_disp, E_elst = calc.get_energy_contributions()
    
    #Pretty print them out
    print(' E_intra = %.4f\n'   % E_intra,
           'E_pair  = %.4f\n'   % E_pair,
           'Sum     = %.4f\n\n' % sum([E_intra, E_pair]))
    print(' E_ex    = %.4f\n'   % E_ex,
           'E_disp  = %.4f\n'   % E_disp,
           'E_elst  = %.4f\n'   % E_elst,
           'Sum     = %.4f\n\n' % sum([E_ex, E_disp, E_elst]))

    #Get analytical forces and track time
    b4     = time.time()
    forces = system.get_forces()
    F_time = time.time() - b4

    #Pretty print analytical forces
    print(' Analytical forces:\n', forces, '\n',
           'Runtime          : %.1f (s)\n\n' % F_time)

    #Get numerical forces and track time
    b4         = time.time()
    num_forces = calc.calculate_numerical_forces(system, d=1e-6)
    F_num_time = time.time() - b4

    #Pretty print numerical forces
    print(' Numerical forces:\n', num_forces, '\n',
           'Runtime         : %.1f (s)\n\n' % F_num_time)

    #Get differences between numerical and analytical
    diff = num_forces - forces
    norm = np.linalg.norm(diff) # also get norm

    #Pretty print difference
    print(' Numerical - Analytical:\n', diff, '\n',
           'Norm                  : %.4f\n\n' % norm )

    #Pretty Closing
    print('\n\n---Ending Printout---\n\n')

    return


