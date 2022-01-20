import sys, time
sys.path.append('../sourcecode/')
from ase.md.verlet import VelocityVerlet
from ase.io import read
from ase import units
from ase import Atoms
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase.vibrations import Vibrations
from MvH_CO_JM8 import MvH_CO
	
#system = Trajectory('moldyn_Large_geo-opt_10K_Langevin_50ps.traj')[-1]
#system = read('Large_isotope_excited.xyz')

#Currently a randomly selected xyz file
system = read('../xyz_files/126.xyz')

#I think here the initial print is to verify the isotope is added
print('\n\nSystem prior to isotope addition:\n')
print(system)

#Add two new atoms (1 CO molecule) but with isotope masses 
#note, it is added on the position of atoms 92 and 93
dt = Atoms("CO", 
           positions=[system.get_positions()[92],system.get_positions()[93]], 
           masses=[13.003, 17.999])

#Delete atoms 92 and 93, then add new isotope atoms 
del system[[92,93]]
system = dt + system
	
#I think this print is to verify the isotope masses are in the system 
print('\n\nSystem after isotope addition:\n')
print(system)

#Add the van Hemet calculator to the system
calc = MvH_CO(atoms=system)
system.set_calculator(calc)
	
#print(system.get_potential_energy())
	
#Initiate the MD simulation with the Verlet numerical method
#Each interval is 1 fs 
dyn = VelocityVerlet(system, 1 * units.fs, 
        logfile='md_Large_isotope_excited_Verlet_10K_50ps.log')

#dyn = Langevin(system, 1*units.fs, 20*units.kB, 0.002, logfile='md_Crystal_geo-opt_20K_Langevin_50ps.log')

 
#Prints the energy/system size
def printenergy(a=system):
    epot = a.get_potential_energy() / len(a)
    print('Potential Energy: ' + str(epot))
	
#Attach a trajectory to the MD, saving every interval
traj = Trajectory('moldyn_Large_isotope_excited_Verlet_10K_50ps.traj', 
                  'w', system)
dyn.attach(traj.write, interval=1)

#Attach the print energy function to the MD, printing ever 100 intervals
dyn.attach(printenergy, interval=100)

#Print the system energy prior to running
printenergy()

#Run for 50k intervals, 50ps 
dyn.run(50000)


#vib = Vibrations(system, delta=0.0001)
#vib.run()
#vib.summary()
#vib.write_dos(out='vib-dos_largeclus_False_w1.dat', start=2175, end=2225, npts=None, width=1, type='Gaussian', method='standard', direction='central')
#vib.write_mode(n=359, kT=0.02585199101165164, nimages=30)


#from ase.optimize import BFGS
#dyn = BFGS(system, trajectory='test.opt.traj')
#dyn.run(fmax=0.0001)
#print(system.get_potential_energy())



