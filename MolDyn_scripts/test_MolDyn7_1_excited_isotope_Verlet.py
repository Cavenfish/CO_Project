#!/usr/bin/env python3

if __name__ == "__main__" :

	import time
	from ase.md.verlet import VelocityVerlet
	from ase.io import read
	from ase import units
	from ase import Atoms
	from ase.md.langevin import Langevin
	from ase.io.trajectory import Trajectory
	
	#system = Trajectory('moldyn_Large_geo-opt_10K_Langevin_50ps.traj')[-1]
	system = read('Large_isotope_excited.xyz')
	from ase.vibrations import Vibrations

	print(system)

	dt = Atoms("CO", positions=[system.get_positions()[92],system.get_positions()[93]], masses=[13.003, 17.999])
	del system[[92,93]]
	system = dt + system
	
	print(system)

	from MvH_CO_JM8 import MvH_CO
	calc = MvH_CO(atoms=system)
	system.set_calculator(calc)
	
	#print(system.get_potential_energy())
	
	
	dyn = VelocityVerlet(system, 1 * units.fs, logfile='md_Large_isotope_excited_Verlet_10K_50ps.log')
	#dyn = Langevin(system, 1*units.fs, 20*units.kB, 0.002, logfile='md_Crystal_geo-opt_20K_Langevin_50ps.log')
	
	def printenergy(a=system):
		epot = a.get_potential_energy() / len(a)
		print(epot)
	
	traj = Trajectory('moldyn_Large_isotope_excited_Verlet_10K_50ps.traj', 'w', system)
	dyn.attach(traj.write, interval=1)
	
	dyn.attach(printenergy, interval=100)
	printenergy()
	dyn.run(50000)
	#vib = Vibrations(system, delta=0.0001)
	#vib.run()
	#vib.summary()
	#vib.write_dos(out='vib-dos_largeclus_False_w1.dat', start=2175, end=2225, npts=None, width=1, type='Gaussian', method='standard', direction='central')
	#vib.write_mode(n=359, kT=0.02585199101165164, nimages=30)


#	from ase.optimize import BFGS
#	dyn = BFGS(system, trajectory='test.opt.traj')
#	dyn.run(fmax=0.0001)
#	print(system.get_potential_energy())



