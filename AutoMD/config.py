import sys, time, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import c
from ase.io import read, write
from ase import units, Atoms
from ase.optimize import BFGS, FIRE, GPMin, MDMin
from .MvH_CO_JM8_Numba import MvH_CO
from ase.visualize import view
from ase.vibrations import Vibrations
from ase.md.verlet import VelocityVerlet
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


def view_xyz(xyz):
    system = read(xyz)
    view(system)
    return

def prep_system(xyz):
    system = read(xyz)
    calc   = MvH_CO(atoms=system)
    system.set_calculator(calc)
    return system, calc

def CoM(pos, masses):
    if not isinstance(pos, np.ndarray):
        pos = np.array(pos)

    M   = sum(masses)
    xcm = sum(masses * pos[:,0])/M
    ycm = sum(masses * pos[:,1])/M
    zcm = sum(masses * pos[:,2])/M

    return np.array([xcm,ycm,zcm])

def rotat_excite(xyz, swap, E):
    #Get pos from system
    system, _ = prep_system(xyz)
    pos       = system.get_positions()[swap]
    momenta   = system.arrays['momenta'][swap]
    masses    = system.get_masses()[swap]

    #Get variables for math
    com = CoM(pos, masses)
    mc  = masses[0]
    mo  = masses[1]
    tc  = pos[0] - com
    rc  = np.sqrt(np.dot(tc,tc))
    to  = pos[1] - com
    ro  = np.sqrt(np.dot(to,to))

    #Get wo and wc
    a   = (mo**2 * ro**2) / mc
    wo  = np.sqrt( (2 * E) / (a + mo * ro**2) )
    wc  = - (mo * wo * ro) / (mc * rc)

    #Get random unitary vector and use it to make velocities
    e  = np.random.rand(3)
    c  = e / np.linalg.norm(e)
    vc = np.cross(wc*c, tc)
    vo = np.cross(wo*c, to)

    #Update the particles momentum
    momenta[0] += mc*vc
    momenta[1] -= mo*vo

    #Make excited molecule
    new_atoms = Atoms('CO', positions=pos, masses=masses, momenta=momenta)

    #Delete old atoms add new atoms
    del system[swap]
    system = new_atoms + system

    #Write XYZ file of system with excited molecule
    new_name = xyz.replace('.xyz', '_excited.xyz')
    write(new_name, system)

    return new_name

def trans_excite(xyz, swap, E):
    #Get pos from system
    system, _ = prep_system(xyz)
    pos       = system.get_positions()[swap]
    momenta   = system.arrays['momenta'][swap]
    masses    = system.get_masses()[swap]

    #Get scalar value of velocity from E
    mu = (masses[0] * masses[1])/sum(masses)
    nu = np.sqrt( (2 * E)/(mu) )

    #Get random unit vector
    r = np.random.rand(3)
    e = r / np.linalg.norm(r)

    #Make velocity vector
    v = nu * e

    #Make new momenta
    p0 = momenta[0] + masses[0] * v
    p1 = momenta[1] + masses[1] * v

    #Make excited molecule
    new_atoms = Atoms('CO', positions=pos, masses=masses, momenta=[p0,p1])

    #Delete old atoms add new atoms
    del system[swap]
    system = new_atoms + system

    #Write XYZ file of system with excited molecule
    new_name = xyz.replace('.xyz', '_excited.xyz')
    write(new_name, system)

    return new_name

def Morse_energy(nu, n):
    #Constants and unit conversions
    D_e   = 11.2301  #eV
    nu    = nu / 1e4 #um^-1
    hc    = 1.239841 #eV um
    tmp   = (n + 1/2) * hc * nu
    V_mor = tmp - ( (tmp**2) / (4 * D_e) )
    return V_mor

def excite_molecule(xyz, swap, E, many_excite=''):
    #Get pos from system
    system, _ = prep_system(xyz)
    pos       = system.get_positions()[swap]
    momenta   = system.arrays['momenta'][swap]
    masses    = system.get_masses()[swap]

    #Define parameters
    R_vect = pos[0] - pos[1]
    R_norm = np.linalg.norm(R_vect)
    R_hat  = R_vect / R_norm

    #Vector points from mol1 to mol0

    #Get new momenta
    a  = 2 * masses[0] * E
    b  = 1 + (masses[0] / masses[1])
    p0 = (np.sqrt(a/b) * R_hat) + momenta[0]
    p1 = momenta[1] - (np.sqrt(a/b) * R_hat)

    #Make excited molecule
    new_atoms = Atoms('CO', positions=pos, masses=masses, momenta=[p0,p1])

    #Delete old atoms add new atoms
    del system[swap]
    system = new_atoms + system

    #Special case for many excite
    if many_excite:
        write(xyz, system)
        return 

    #Write XYZ file of system with isotope
    new_name = xyz.replace('.xyz', '_excited.xyz')
    write(new_name, system)

    return new_name

def excite_many_molecules(xyz, n, E):
    system, _ = prep_system(xyz)
    N         = len(system)
    new_name  = xyz.replace('.xyz', '_' + str(n) + 'excited.xyz')
    write(new_name, system)

    for i in range(n+1):
        j = np.random.randint(0, N//2)
        excite_molecule(new_name, [j*2, j*2+1], E, new_name)

    return new_name

def stretch_molecule(xyz, swap, r):
    #r excitation radius

    #Get pos from system
    system, _ = prep_system(xyz)
    pos       = system.get_positions()[swap]
    momenta   = system.arrays['momenta'][swap]
    masses    = system.get_masses()[swap]

    #Define parameters
    R_vect = pos[0] - pos[1]
    R_norm = np.linalg.norm(R_vect)
    R_prim = r * (R_vect / R_norm)

    top    = masses[0]*pos[0] + masses[1]*pos[1] + masses[0]*R_prim
    bottom = sum(masses)
    r2p    = top/bottom
    r1p    = R_prim + r2p

    new_pos  = [r1p, r2p]

    #Make stretched molecule
    new_atoms = Atoms('CO', positions=new_pos, masses=masses, momenta=momenta)

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
    D_e   = 11.2301  #eV
    beta  =  2.627   #Angstrom^-1
    r_e   =  1.1282  #Angstrom
    nu    = nu / 1e4   #um^-1
    hc    = 1.239841   #eV um
    tmp   = (n + 1/2) * hc * nu
    V_mor = tmp - ( (tmp**2) / (4 * D_e) )

    #V_mor = D_e[1-exp(-beta(r_A-r_e))]^2

    a = np.sqrt(V_mor/D_e)

    #a = 1-exp(-beta(r_A-r_e)) ---> exp(-beta(r_A-r_e)) = 1 - a

    b = - np.log(1-a)/beta

    #r_A - r_e = b ----> r_A = b + r_e

    r_A = b + r_e

    return r_A

def geo_opt(xyz, fmax=0.0001, method='BFGS'):
    #Read in system and set van Hemert calculator
    system, calc = prep_system(xyz)

    #Make trajectory string
    traj  = xyz.replace('.xyz', '_opt.traj')

    #Methods dictionary
    meth  = {'BFGS' :  BFGS(system, trajectory=traj),
             'FIRE' :  FIRE(system, trajectory=traj),
             'GPMin': GPMin(system, trajectory=traj),
             'MDMin': MDMin(system, trajectory=traj)}

    #Run optimization of geometry
    opt   = meth[method]
    opt.run(fmax=fmax)

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
    vib.summary(log=xyz.replace('.xyz', '_vib.out'))
    vib.write_dos(xyz.replace('.xyz', '_vibDOS.out'))
    vib.write_jmol()
    vib.clean()
    return

def add_isotope(xyz, swap, masses):
    #Read in system
    system, _ = prep_system(xyz)

    #Get atom positions
    pos       = system.get_positions()[swap]
    momenta   = system.arrays['momenta'][swap]

    #Make new atoms but with isotopic masses
    new_atoms = Atoms('CO', positions=pos, masses=masses, momenta=momenta)

    #Delete old atoms add new atoms
    del system[swap]
    system = new_atoms + system

    #Write XYZ file of system with isotope
    new_name = xyz.replace('.xyz', '_isotope.xyz')
    write(new_name, system)

    return new_name

def run_langevinMD(xyz, n=50000, mu=0.002, temp=10, i=1000):
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
    dyn.attach(traj.write, interval=i)

    #Run for n intervals (1 fs/interval)
    dyn.run(n)

    #Get last frame of simulation write xyz for it
    traj  = Trajectory(trajname)
    atoms = traj[-1]
    fname = xyz.replace('.xyz', '_NVT.xyz')
    write(fname, atoms)

    return fname

def run_verletMD(xyz, n=50000, i=100, ts=1):
    #Prep system
    system, calc = prep_system(xyz)

    #Define logfile name
    logfile  = xyz.replace('.xyz', '_NVE.log')

    #Initiate MD simulation with Verlet numerical method
    dyn      = VelocityVerlet(system, ts * units.fs, logfile=logfile)

    #Attach a trajectory file to the MD, saving every interval
    trajname = xyz.replace('.xyz', '_NVE.traj')
    traj     = Trajectory(trajname, 'w', system)
    dyn.attach(traj.write, interval=i)

    #Run for n intervals (1 fs/interval)
    dyn.run(n)

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

def track_dissipation(system, calc):
    #Get potential energy of molecule
    try:
        E_pot   = system.calc.results['energies'][0]
    except:
        E_pot   = 0

    #Get kinetic energy of molecule
    momenta = system.arrays['momenta']
    masses  = system.arrays['masses' ]
    E_kin_a = (np.dot(momenta[0], momenta[0])) / (2 * masses[0])
    E_kin_b = (np.dot(momenta[1], momenta[1])) / (2 * masses[1])
    E_kin   = E_kin_a + E_kin_b

    #Sum up total molecular energy
    E_tot   = E_kin + E_pot

    #Get total energy for sliced system
    part  = system[0:2]
    part.set_calculator(MvH_CO(atoms=part))
    E_sli = part.get_potential_energy() + part.get_kinetic_energy()

    #Get kinetic energy contributions of excited molecule
    pos      = system.arrays['positions']
    velocs   = system.get_velocities()
    one_vib  = calc_Evib(pos[0:2], masses[0:2], velocs[0:2])
    one_rot  = calc_Erot(pos[0:2], masses[0:2], velocs[0:2])
    one_tran = calc_Etran(pos[0:2], masses[0:2], velocs[0:2])

    #Get kinetic energy contributions of all other molecules
    N = len(pos) // 2
    all_vib  = []
    all_rot  = []
    all_tran = []
    for i in range(1, N):
        a = i*2
        b = a+2
        all_vib.append(calc_Evib(pos[a:b], masses[a:b], velocs[a:b]))
        all_rot.append(calc_Erot(pos[a:b], masses[a:b], velocs[a:b]))
        all_tran.append(calc_Etran(pos[a:b], masses[a:b], velocs[a:b]))
    avg_vib  = sum(all_vib)  / (N-1)
    avg_rot  = sum(all_rot)  / (N-1)
    avg_tran = sum(all_tran) / (N-1)

    return_this = [E_sli, E_pot, E_kin, E_tot,
                   avg_tran, avg_rot, avg_vib,
                   one_tran, one_rot, one_vib,
                   all_tran, all_rot, all_vib]
    return return_this

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

def make_NVE_output(trajFile, csvFile, ts, mask=False):
    #Open trajectory
    traj = Trajectory(trajFile)

    #Initiate dictionary
    temp = {            'Time': [],
               'Sliced Energy': [],
            'Potential Energy': [],
              'Kinetic Energy': [],
                'Total Energy': [],
            'Avg Trans Energy': [],
            'Avg Rotat Energy': [],
            'Avg Vibra Energy': [],
            'One Trans Energy': [],
            'One Rotat Energy': [],
            'One Vibra Energy': [],
            'All Trans Energy': [],
            'All Rotat Energy': [],
            'All Vibra Energy': []}

    #Loop through trajectory, writting to output file
    for i in range(len(traj)):
        system = traj[i]
        if mask:
            rList  = np.array(get_rList(system))
            system = system[(mask[0] <= rList) & (rList <= mask[1])]
        calc   = system.calc
        data   = track_dissipation(system, calc)
        temp['Time'].append(i * ts)
        temp['Sliced Energy'].append(data[0])
        temp['Potential Energy'].append(data[1])
        temp['Kinetic Energy'].append(data[2])
        temp['Total Energy'].append(data[3])
        temp['Avg Trans Energy'].append(data[4])
        temp['Avg Rotat Energy'].append(data[5])
        temp['Avg Vibra Energy'].append(data[6])
        temp['One Trans Energy'].append(data[7])
        temp['One Rotat Energy'].append(data[8])
        temp['One Vibra Energy'].append(data[9])
        temp['All Trans Energy'].append(data[10])
        temp['All Rotat Energy'].append(data[11])
        temp['All Vibra Energy'].append(data[12])

    #Make DataFrame
    df = pd.DataFrame(temp)

    #Write csv file
    df.to_csv(csvFile)

    #Return DataFrame in case its going to be used
    #If not, the call can ignore the return
    return df

def calc_Evib(pos, masses, velocs):
    #Calculate kinetic energy terms
    d     = pos[0] - pos[1]
    ed    = d / np.linalg.norm(d)
    vd0   = np.dot(ed, velocs[0])
    vd1   = np.dot(ed, velocs[1])
    a     = 0.5 * masses[0] * vd0**2
    b     = 0.5 * masses[1] * vd1**2

    #Sum it all up
    E_vib = a + b
    return E_vib

def calc_Erot(pos, masses, velocs):
    vCoM  = (masses[0]*velocs[0] + masses[1]*velocs[1]) / sum(masses)
    com   = CoM(pos, masses)
    E_rot = 0
    for i in range(len(masses)):
        r = pos[i] - com
        v = velocs[i] - vCoM
        m = masses[i]
        w = np.cross(r, v) / np.dot(r, r)
        I = m * np.dot(r, r)

        E_rot += 0.5 * I * np.dot(w, w)
    return E_rot

def calc_Etran(pos, masses, velocs):
    mu     = (masses[0] * masses[1]) / sum(masses)
    vCoM   = (masses[0]*velocs[0] + masses[1]*velocs[1]) / sum(masses)
    E_tran = 0.5 * mu * np.dot(vCoM, vCoM)
    return E_tran

def get_spherical_coords(xyz):
    system = read(xyz)
    pos    = system.arrays['positions']
    masses = system.arrays['masses']
    com    = CoM(pos, masses)
    coords = []

    for i in range(len(masses) // 2):
        a     = i*2
        b     = a+2
        cm    = CoM(pos[a:b], masses[a,b])
        x,y,z = cm
        diff  = cm - com
        r     = np.sqrt(np.dot(diff, diff))
        theta = np.arccos(z/r)
        phi   = np.arctan2(y,x)

        coords.append([r, theta, phi])

    return coords

def get_rList(atoms):
    pos    = atoms.get_positions()
    masses = atoms.get_masses()
    eCoM   = CoM(pos[:2], masses[:2])
    rList  = []
    for i in range(0, len(masses)//2):
        a     = i*2
        b     = a+2
        cm    = CoM(pos[a:b], masses[a:b])
        diff  = cm - eCoM
        r     = np.sqrt(np.dot(diff,diff))
        rList.append(r)
        rList.append(r)

    return rList

def track_E_vs_r(EList, pos, masses):
    eCoM   = CoM(pos[:2], masses[:2])

    rList  = []
    for i in range(1, len(masses)//2):
        a     = i*2
        b     = a+2
        cm    = CoM(pos[a:b], masses[a:b])
        diff  = cm - eCoM
        r     = np.sqrt(np.dot(diff,diff))
        rList.append(r)

    plt.scatter(rList, EList)
    #plt.show()

    return rList

def radial_distribution(xyz):
    atoms = read(xyz)

    pos    = atoms.get_positions()
    masses = atoms.get_masses()
    eCoM   = CoM(pos[:2], masses[:2])
    C      = pos[0]
    O      = pos[1]
    coco   = []
    co     = []
    cc     = []
    oo     = []
    for i in range(1, len(masses)//2):
        a     = i*2
        b     = a+2
        cm    = CoM(pos[a:b], masses[a:b])
        diff  = cm - eCoM
        r     = np.sqrt(np.dot(diff,diff))
        coco.append(r)
        c     = pos[a]
        o     = pos[a+1]
        r     = np.sqrt(np.dot(C-o, C-o))
        co.append(r)
        r     = np.sqrt(np.dot(C-c, C-c))
        cc.append(r)
        r     = np.sqrt(np.dot(O-o, O-o))
        oo.append(r)

    return coco,co,cc,oo
