
import numpy as np
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.calculators.calculator import Calculator
traj = Trajectory('126_FQ_False_JM7_fmax_e-4.opt.traj')
atoms = traj[-1]
print(type(atoms))
print(traj)
print(atoms)
positions = atoms.get_positions()
print(positions)

#with open("final_config_cluster2_FQ_MODE_True_fmaxE-4.txt", "w") as f:
#	print(positions, file=f)
N = len(positions)
for i in range(0, N//2, 1):
	difference = positions[2*i+1]-positions[2*i]
	distance = np.linalg.norm(difference)
	print(difference, distance,i)
	with open("CO_distances_CO_crystal_False_fmaxE-4.txt", 'a') as f:
		print(distance, file=f)




	

