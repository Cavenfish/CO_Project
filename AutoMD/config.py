import sys, time
import numpy as np
from scipy.constants import c
from ase.io import read, write
from ase import units, Atoms
from ase.optimize import BFGS
from .MvH_CO_JM8 import MvH_CO
from ase.visualize import view
from ase.vibrations import Vibrations
from ase.md.verlet import VelocityVerlet
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin

def view_xyz(xyz):
    system = read(xyz)
    view(system)
    return 

def prep_system(xyz):
    system = read(xyz)
    calc   = MvH_CO(atoms=system)
    system.set_calculator(calc)
    return system, calc 

def stretch_molecule(xyz, swap, masses, pos, r):
    #r excitation radius

    #concept: reduce the problem to 1D

    #CoM = (m_c*x_c + m_o*x_o)/(m_c + m_o)
    # l  = |x_c - x_o|

    #idea: -Find vector from C->O
    #      -Make unitary vector of it (this will act as hat{x})
    #      -Move C and O along opposite directions of pseudo hat{x}
    #      -Weight the movement based on their mass ratio

    bnd_vect = np.array(pos[1]) - np.array(pos[0])
    norm     = np.linalg.norm(bnd_vect)
    unt_vect = bnd_vect / norm
    r_diff   = r - norm
    move     = unt_vect*r_diff

    #note: unit vector points from pos[0] to pos[1]
    #ie. positive goes toward pos[1]
    #    negative goes away from pos[1]
    
    #note: I weight the movement by the other atoms 
    #      mass ratio since the displacement is larger
    #      for smaller weight rather than larger weight
    
    m_ratio0 = masses[1] / sum(masses)
    m_ratio1 = masses[0] / sum(masses)
    new_pos0 = np.array(pos[0]) - (move * m_ratio0)
    new_pos1 = np.array(pos[1]) + (move * m_ratio1)
    
    new_pos  = [new_pos0, new_pos1]

    #Read in xyz file
    system = read(xyz)
    
    #Make stretched molecule
    new_atoms = Atoms('CO', positions=new_pos, masses=masses)

    #Delete old atoms add new atoms
    del system[swap]
    system = new_atoms + system

    #Write XYZ file of system with isotope
    new_name = xyz.replace('.xyz', '_excited.xyz')
    write(new_name, system)

    return

def Morse_excitation(nu, n):
    #nu -> cm^-1
    # n -> unitless

    #Note: D_e, r_e and beta are values fitted to CO
    #      if using an isotope, we need new values

    #Constants and unite conversions
    D_e   = 11.2301       #eV
    beta  =  0.6519       #Angstrom^-1
    r_e   =  1.1282       #Angstrom
    omega = nu * c * 100  #s^-1
    hbar  =  6.582119e-16 #eV * s
    V_mor = (n + 1/2) * hbar * omega

    #V_mor = D_e[1-exp(-beta(r_A-r_e))]^2

    a = np.sqrt(V_mor/D_e)

    #a = 1-exp(-beta(r_A-r_e)) ---> exp(-beta(r_A-r_e)) = 1 - a

    b = - np.log(1-a)/beta

    #r_A - r_e = b ----> r_A = b + r_e

    r_A = b + r_e

    return r_A

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
    positions = [system.get_positions()[pos[0]], 
                 system.get_positions()[pos[1]]]

    #Make new atoms but with isotopic masses
    new_atoms = Atoms('CO', positions=positions, masses=masses)

    #Delete old atoms add new atoms
    del system[pos]
    system = new_atoms + system

    #Write XYZ file of system with isotope
    new_name = xyz.replace('.xyz', '_isotope.xyz')
    write(new_name, system)

    return 

def run_langevinMD(xyz, n=50000, pos=False, masses=False):
    #If posisitons and masses given, add isotope and prep system
    #Else just prep given system
    if (pos) and (masses):
        system, calc = add_isotope(xyz, pos, masses)
        xyz          = xyz.replace('.xyz', '_isotope.xyz')
    else:
        system, calc = prep_system(xyz)

    #Define logfile name
    logfile  = xyz.replace('.xyz', '.log')

    #Initiate MD simulation with Langevin thermometer
    dyn      = Langevin(system, 
                        1  * units.fs, #time interval
                        20 * units.kB, #temperature
                        0.002,         #friction coefficient
                        logfile=logfile)

    #Attach a trajectory file to the MD, saving every interval
    trajname = xyz.replace('.xyz', '.traj')
    traj     = Trajectory(trajname, 'w', system)
    dyn.attach(traj.write, interval=1)

    #Open output file
    outname = xyz.replace('.xyz', '.out')
    out     = open(outname, 'w')
    f = lambda x=system, y=out: (
            y.write(str(x.get_potential_energy() / len(x)) +'\n'))

    #Attach the lambda function to the MD, every 100 intervals
    dyn.attach(f, interval=100)

    #Run for n intervals (1 fs/interval)
    dyn.run(n)

    #Save and close output file
    out.close()

    return

def run_verletMD(xyz, n=50000,  pos=False, masses=False):
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

    #Open output file
    outname = xyz.replace('.xyz', '.out')
    out     = open(outname, 'w')
    f = lambda x=system, y=out: (
            y.write(str(x.get_potential_energy() / len(x)) +'\n'))

    #Attach the lambda function to the MD, every 100 intervals
    dyn.attach(f, interval=100)

    #Run for n intervals (1 fs/interval)
    dyn.run(n)

    #Save and close output file
    out.close()

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


