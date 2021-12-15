#!/usr/bin/env python3

if __name__ == "__main__" :

	import time
	from ase.io.trajectory import Trajectory
	from ase.io import read
	system = read('optimized_pentamer.xyz')

	from ase.vibrations import Vibrations

	from MvH_CO_JM8 import MvH_CO
#	from coco_MvH_ASE import coco_MvH
#	from coco_MvH_ASE import coco_qfix_MvH
#	from coco_MvH_ASE import coco2_MvH
	calc = MvH_CO(atoms=system)
	system.set_calculator(calc)


	
	print(system.get_potential_energy())
	if isinstance(calc, MvH_CO) or isinstance(calc, coco2_MvH):
		E_intra, E_pair, E_ex, E_disp, E_elst = calc.get_energy_contributions()
		print(sum([E_intra,E_pair]), E_intra, E_pair)
		print(sum([E_ex, E_disp, E_elst]), E_ex, E_disp, E_elst)

	start = time.time()
	analytical_forces = system.get_forces()
	end = time.time()
	print("Analytical forces: ", end-start, "seconds \n", analytical_forces)

	start = time.time()
	numerical_forces = calc.calculate_numerical_forces(system, d=1E-6)
	end = time.time()
	print("Numerical forces: ", end-start, "seconds \n", numerical_forces)
	
	difference = numerical_forces - analytical_forces
	print("Difference:")
	print(difference)
	import numpy as np
	verschil = np.linalg.norm(numerical_forces - analytical_forces)
	print(verschil)

	
	vib = Vibrations(system, delta=0.0001)
	vib.run()
	vib.summary()


#	from ase.optimize import BFGS
#	dyn = BFGS(system, trajectory='Pentamer_FQ_True_JM7_fmax_e-4.opt.traj')
#	dyn.run(fmax=0.0001)
#	print(system.get_potential_energy())

#	with open("Pentamer_FQ_True_JM7_fmax_e-4.txt", "a") as f:
	
#		print(system.get_potential_energy(), file=f)
#		if isinstance(calc, MvH_CO) or isinstance(calc, coco2_MvH):
#			E_intra, E_pair, E_ex, E_disp, E_elst = calc.get_energy_contributions()
#			print(sum([E_intra,E_pair]), E_intra, E_pair, file=f)
#			print(sum([E_ex, E_disp, E_elst]), E_ex, E_disp, E_elst, file=f)

#		start = time.time()
#		analytical_forces = system.get_forces()
#		end = time.time()
#		print("Analytical forces: ", end-start, "seconds \n", analytical_forces, file=f)

#		start = time.time()
#		numerical_forces = calc.calculate_numerical_forces(system, d=1E-6)
#		end = time.time()
#		print("Numerical forces: ", end-start, "seconds \n", numerical_forces, file=f)
	
#		difference = numerical_forces - analytical_forces
#		print("Difference:", file=f)
#		print(difference, file=f)
#		import numpy as np
#		verschil = np.linalg.norm(numerical_forces - analytical_forces)
#		print(verschil, file=f)
	# cross-check energy of local minimum with a differen calculator
# 	calc = coco2_MvH(atoms=system)
# 	system.set_calculator(calc)
# 	print(system.get_potential_energy())
