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

def add_isotope(xyz, pos, masses):
    #Read in system
    system    = read(xyz)

    #Get atom positions
    positions = [system.get_positions()[pos]]

    #Make new atoms but with isotopic masses
    new_atoms = Atoms('CO', positions=positions, masses=masses)

    #Delete old atoms add new atoms
    del system[pos]
    system = new_atoms + system

    #Write XYZ file of system with isotope
    new_name = xyz.replace('.xyz', '_isotope.xyz')
    write(new_name, system)

    return system, calc 


def run_verletMD(xyz, pos=False, masses=False):
    #If posisitons and masses given, add isotope and prep system
    #Else just prep given system
    if (pos) and (masses):
        system, calc = add_isotope(xyz, pos, masses)
        xyz          = xyz.replace('.xyz', '_isotope.xyz')
    else:
        system, calc = prep_system(xyz)

    #Define logfile name
    logfile  = xyz.replace('.xyz', '.log') 

    #Initiate MD simulation with Verlet numerical method
    dyn      = VelocityVerlet(system, 1 * units.fs, logfile=logfile)

    #Attach a trajectory file to the MD, saving every interval
    trajname = xyz.replace('.xyz', '.traj')
    traj     = Trajectory(trajname, 'w', system)
    dyn.attach(traj.write, interval=1)

    f = lambda x=system: (print(x.get_potential_energy() / len(x))) 

    #Attach the lambda function to the MD, every 100 intervals
    dyn.attach(f, interval=100)

    #Run for 50k intervals (1 fs/interval -> 50 ps total)
    dyn.run(50000)

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


