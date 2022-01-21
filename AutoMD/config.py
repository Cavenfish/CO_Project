import sys, time
from ase.io import read
from ase import units, Atoms
from ase.optimize import BFGS
from MvH_CO_JM8 import MvH_CO
from ase.visualize import view
from ase.vibrations import Vibrations
from ase.md.verlet import VelocityVerlet

def view_xyz(xyz):
    system = read(xyz)
    view(system)
    return 

def prep_system(xyz):
    system = read(xyz)
    calc   = MvH_CO(atoms=system)
    system.set_calculator(calc)
    return system

def geo_opt(xyz):
    #Read in system and set van Hemet calculator
    system = prep_system(xyz)

    #Make trajectory string
    traj   = xyz.replace('.xyz', '_opt.traj')
    
    #Run BFGS optimization of geometry
    opt    = BFGS(system, trajectory=traj)
    opt.run(fmax=0.0001)

    return

def get_system_properties(xyz):
    #Read in system and set van Hemet calculator
    system = prep_system(xyz)

    #Get system potential energy and track time of calculation
    b4          = time.time()
    E_pot       = system.get_potential_energy())
    E_pot_time  = time.time() - b4
    set0        = (E_pot, E_pot_time)

    #Pretty print potential energy
    print('Potential Energy: %.4f | Runtime: %.1f (s)\n\n' % set0)

    #Get energy contributions
    E_intra, E_pair, E_ex, E_disp, E_elst = calc.get_energy_contributions()
    
    #Pretty print them out
    set1 = (E_intra, E_pair, sum([E_intra,E_pair]))
    set2 = (E_ex, E_disp, E_elst, sum([E_ex, E_disp, E_elst]))
    print(' E_intra = %.4f\n E_pair  = %.4f\n Sum     = %.4f\n\n' % set1)
    print(' E_ex    = %.4f\n E_disp  = %.4f\n E_elst  = %.4f\n' +
          ' Sum     = %.4f\n\n' %set2)

    #Get forces and track time
    b4     = time.time()
    forces = system.get_forces()
    F_time = time.time() - b4
    set3   = (forces, F_time)

    #Pretty print forces
    print('Analytical forces: %.4f | Runtime: %.1f (s)\n\n' % set3)

    return
