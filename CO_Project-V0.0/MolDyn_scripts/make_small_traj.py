from ase.io import read
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.geometry.analysis import Analysis
from ase import Atoms

system1 = Trajectory('moldyn_Large_excited_Verlet_20Ksteps_10K_20ps.traj')
system2 = Trajectory('moldyn_Large_excited_Verlet_30Ksteps_10K_50ps.traj')
atoms_list = []

for i in range(len(system1)):
	CO_positions = system1[i][92:94].get_positions()
	CO_velocities = system1[i][92:94].get_velocities()
	CO_append = Atoms(symbols='CO', positions=CO_positions, velocities = CO_velocities)
#	print(CO_append)
	atoms_list.append(CO_append)

for i in range(len(system2)):
	CO_positions = system2[i][92:94].get_positions()
	CO_velocities = system2[i][92:94].get_velocities()
	CO_append = Atoms(symbols='CO', positions=CO_positions, velocities = CO_velocities)
	atoms_list.append(CO_append)

print(atoms_list)
new_traj = write('CO_molecule_50ps.traj', images = atoms_list)
	
print(new_traj)

#test = Trajectory('CO_molecule.traj')
#print(type(test))

