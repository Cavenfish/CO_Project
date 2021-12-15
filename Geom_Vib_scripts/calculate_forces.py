#!/usr/bin/env python3

if __name__ == "__main__" :

	import time

	from ase.io import read
	system = read('pentamer.xyz')

	from MvH_CO_JM8 import MvH_CO
	calc = MvH_CO(atoms=system)
	system.set_calculator(calc)

#	with open("Large_cluster_FQ_False_JM7_fmax_e-4.txt", "a") as f:
	b4 = time.time()
	print(system.get_potential_energy())#, file=f)
	after_E_pot = time.time()
	if isinstance(calc, MvH_CO) or isinstance(calc, coco2_MvH):
		E_intra, E_pair, E_ex, E_disp, E_elst = calc.get_energy_contributions()
		print(sum([E_intra,E_pair]), E_intra, E_pair)#, file=f)
		print(sum([E_ex, E_disp, E_elst]), E_ex, E_disp, E_elst)#, file=f)

	start = time.time()
	analytical_forces = system.get_forces()
	end = time.time()
	print("Analytical forces: ", end-start, "seconds \n", analytical_forces)#, file=f)

	start = time.time()
#	numerical_forces = calc.calculate_numerical_forces(system, d=1E-6)
#	end = time.time()
#	print("Numerical forces: ", end-start, "seconds \n", numerical_forces)#, file=f)
	
#	difference = numerical_forces - analytical_forces
#	print("Difference:")#, file=f)
#	print(difference)#, file=f)
#	import numpy as np
#	verschil = np.linalg.norm(numerical_forces - analytical_forces)
#	print(verschil)#, file=f)

	time_for_E_pot = after_E_pot - b4
	print('E pot time is:',time_for_E_pot)

	from ase.optimize import BFGS
	dyn = BFGS(system, trajectory='pentamer_opt.traj')
	dyn.run(fmax=0.0001)
	print(system.get_potential_energy())

#	with open("Large_cluster_FQ_False_JM7_fmax_e-4.txt", "a") as f:
	
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
