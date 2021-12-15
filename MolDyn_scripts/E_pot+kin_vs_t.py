from ase.io import read
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.geometry.analysis import Analysis
from ase import Atoms
from MvH_CO_JM8 import MvH_CO
import matplotlib.pyplot as plt


CO_traj = Trajectory('CO_molecule_50ps.traj')
E_pot_list = []
E_kin_list = []
E_total = []
for i in range(len(CO_traj)):
	system = CO_traj[i]	
	calc = MvH_CO(atoms=system)
	system.set_calculator(calc)
	E_pot = system.get_potential_energy()
	E_pot_list.append(E_pot)
	E_kin = system.get_kinetic_energy()
	E_kin_list.append(E_kin)

for i in range(len(E_pot_list)):
	E_total.append(E_pot_list[i] + E_kin_list[i])

E_total_avg = []

for i in range(len(E_total)//200):
	E_total_interval = 0
	for j in range(200):
		E_total_interval += E_total[200*i+j]
	E_total_avg.append(E_total_interval/200)
#print(E_total_avg)

#print(E_total)

#plt.plot(E_pot_list)
#plt.plot(E_kin_list)
plt.plot(range(0,50000,200),E_total_avg)
plt.ylabel('$E_{Tot}$(eV)')
plt.xlabel('Simulation time(fs)')
#print(E_pot_list)
plt.show()
