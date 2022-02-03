import sys, time
import numpy as np
import pandas as pd
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

def stretch_molecule(xyz, swap, masses, r):
    #r excitation radius

    #concept: reduce the problem to 1D

    #CoM = (m_c*x_c + m_o*x_o)/(m_c + m_o)
    # l  = |x_c - x_o|

    #idea: -Find vector from C->O
    #      -Make unitary vector of it (this will act as hat{x})
    #      -Move C and O along opposite directions of pseudo hat{x}
    #      -Weight the movement based on their mass ratio

    #Get pos from system
    atoms, _ = prep_system(xyz)
    pos      = atoms.get_positions()[swap]

    #Define parameters 
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

    return new_name

def Morse_excitation(nu, n):
    #nu -> cm^-1
    # n -> unitless

    #Note: D_e, and beta are values fitted to CO
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

    return opt_f

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

    return new_name

def run_langevinMD(xyz, n=50000, mu=0.002, temp=10):
    #Prep system
    system, calc = prep_system(xyz)

    #Define logfile name
    logfile  = xyz.replace('.xyz', '_NVT.log')

    #Initiate MD simulation with Langevin thermometer
    dyn      = Langevin(system, 
                        1  * units.fs, #time interval
                        friction      = mu,
                        temperature_K = temp,
                        logfile       =logfile)

    #Attach a trajectory file to the MD, saving every interval
    trajname = xyz.replace('.xyz', '_NVT.traj')
    traj     = Trajectory(trajname, 'w', system)
    dyn.attach(traj.write, interval=1)

    #UPDATE: This adds unneccsary time to MD runtime
    #            I am removing it
    #Open output file
    #outname = xyz.replace('.xyz', '_NVT.out')
    #out     = open(outname, 'w')
    #f = lambda x=system, y=out: (
    #        y.write(str(x.get_potential_energy() / len(x)) +'\n'))

    #Attach the lambda function to the MD, every 100 intervals
    #dyn.attach(f, interval=1)

    #Run for n intervals (1 fs/interval)
    dyn.run(n)

    #Save and close output file
    #out.close()

    #Get last frame of simulation write xyz for it
    traj  = Trajectory(trajname)
    atoms = traj[-1]
    fname = xyz.replace('.xyz', '_NVT.xyz')
    write(fname, atoms)

    return fname

def run_verletMD(xyz, n=50000):
    #Prep system
    system, calc = prep_system(xyz)

    #Define logfile name
    logfile  = xyz.replace('.xyz', '_NVE.log') 

    #Initiate MD simulation with Verlet numerical method
    dyn      = VelocityVerlet(system, 1 * units.fs, logfile=logfile)

    #Attach a trajectory file to the MD, saving every interval
    trajname = xyz.replace('.xyz', '_NVE.traj')
    traj     = Trajectory(trajname, 'w', system)
    dyn.attach(traj.write, interval=1)


    #UPDATE: This adds unneccsary time to MD runtime
    #            I am removing it

    #Open output file
    #outname = xyz.replace('.xyz', '_NVE.out')
    #out     = open(outname, 'w')
    #f = lambda x=system: (x.get_potential_energy())

    #Attach the lambda function to the MD, every 100 intervals
    #dyn.attach(f, interval=1)

    #Run for n intervals (1 fs/interval)
    dyn.run(n)

    #Save and close output file
    #out.close()

    return trajname

def get_system_properties(system, calc, i, f):
    #Get system potential energy and track time of calculation
    b4          = time.time()
    E_pot       = system.get_potential_energy()
    #E_pots      = system.get_potential_energies()
    E_pot_time  = time.time() - b4

    #Pretty Header
    f.write('\n\n---Interval %d Printout---\n\n' % i)

    #Pretty print excited molecule energy
    #f.write(' Excited Molecule Energy: %.4f\n\n' % sum(E_pots[0:2]))

    #Pretty print potential energy
    f.write('Potential Energy: %.4f\n' % E_pot)
    f.write('Runtime         : %.1f (s)\n\n' % E_pot_time)

    #Get energy contributions
    E_intra, E_pair, E_ex, E_disp, E_elst = calc.get_energy_contributions()

    #Pretty print them out
    f.write('E_intra = %.4f\n'   % E_intra)
    f.write('E_pair  = %.4f\n'   % E_pair)
    f.write('Sum     = %.4f\n\n' % sum([E_intra, E_pair]))
    f.write('E_ex    = %.4f\n'   % E_ex)
    f.write('E_disp  = %.4f\n'   % E_disp)
    f.write('E_elst  = %.4f\n'   % E_elst)
    f.write('Sum     = %.4f\n\n' % sum([E_ex, E_disp, E_elst]))

    #Get analytical forces and track time
    b4     = time.time()
    forces = system.get_forces()
    F_time = time.time() - b4

    #Pretty print analytical forces
    f.write('Analytical forces:\n' +  str(forces) + '\n')
    f.write('Runtime          : %.1f (s)\n\n' % F_time)

    #Get numerical forces and track time
    b4         = time.time()
    num_forces = calc.calculate_numerical_forces(system, d=1e-6)
    F_num_time = time.time() - b4

    #Pretty print numerical forces
    f.write('Numerical forces:\n' + str(num_forces) + '\n')
    f.write('Runtime         : %.1f (s)\n\n' % F_num_time)

    #Get differences between numerical and analytical
    diff = num_forces - forces
    norm = np.linalg.norm(diff) # also get norm

    #Pretty print difference
    f.write('Numerical - Analytical:\n' + str(diff) + '\n')
    f.write('Norm                  : %.4f\n\n' % norm )

    #Pretty Closing
    f.write('\n\n---Ending Printout---\n\n')

    return

def track_dissipation(system, calc, i, df):
    #Idea for tracking single molecule energy:
    #   -Get energy of full system
    #   -Get energy of system minus molecule (slice system)
    #   -Full system - sliced system =  single molecule + disp
    #   -Single molecule + disp is what i want

    #UPDATE: This method is extremly slow!
    #           it takes roughly 0.07s per frame for 17CO cluster
    #           it takes roughly 0.85s per frame for 60CO cluster

    #TEMP SOLU: Return to using the single molecule slice
    #               this won't include the effect of the
    #               other molecules but should still give
    #               some information. 

    #Fastest method is to store molecule energ during MD
    #Next best is calculate just molecular energy with 
    #specific function that uses pair potential

    #Set up a partial system Atoms object with calculator

    #Get potential energy of molecule
    #E_pot_full  =  system.get_potential_energy()
    #E_pot_rest  = partial.get_potential_energy()
    
    #NEW METHOD: Energy of molecules is saved during MD 
    #               now we just pull that saved data
    #               should be nearly 0 runtime per interval
    #
    #NOTICE:     I am switching to momenta based E_kin formula
    #               this is because velocities need to be calculated
    #               whereas the momenta is pulled from saved space

    #Get potential energy of molecule
    E_pot   = system.calc.results['energies'][0]

    #Get kinetic energy of molecule 
    momenta = system.arrays['momenta']
    masses  = system.arrays['masses' ]
    E_kin_a = (np.dot(momenta[0], momenta[0])) / (2 * masses[0])
    E_kin_b = (np.dot(momenta[0], momenta[0])) / (2 * masses[0])
    E_kin   = E_kin_a + E_kin_b 
    
    #Sum up total molecular energy
    E_tot   = E_kin + E_pot

    #Initiate row of data to populate DataFrame with
    new_row = {'Time': i,
               'Potential Energy': E_pot,
               'Kinetic Energy': E_kin,
               'Total Energy': E_tot}
    
    #Add data to DataFrame
    df = df.append(new_row, ignore_index=True)

    return df

def make_NVT_output(logFile, csvFile):
    #Open, read and close log file
    f    = open(logFile, 'r')
    data = f.readlines()
    f.close()

    #Get header row and its length
    head = data[0].split()
    
    #Initiate dictionary
    temp = {}

    #Nest loops to build dict
    for i in range(len(head)):
        temp[head[i]] = []
        for j in range(1,len(data)):
            temp[head[i]].append(float(data[j].split()[i]))

    #Build DataFrame 
    df = pd.DataFrame(temp)

    #Make csv file
    df.to_csv(csvFile)

    #Return DataFrame in case its going to be used
    #If not, the call can ignore the return
    return df 

def make_NVE_output(trajFile, csvFile):
    #Open trajectory
    traj = Trajectory(trajFile)
    
    #Initiate DataFrame
    df   = pd.DataFrame({'Time':[],
                         'Potential Energy': [],
                         'Kinetic Energy': [],
                         'Total Energy': []})

    #NOTICE: Current runtime per loop is 0.002s
    #           might be able to speed it up by 
    #           filling a dict first then making a DF

    #Loop through trajectory, writting to output file
    for i in range(len(traj)):
        system = traj[i]
        calc   = system.calc
        df     = track_dissipation(system, calc, i, df)

    #Write csv file
    df.to_csv(csvFile)

    #Return DataFrame in case its going to be used
    #If not, the call can ignore the return
    return df
