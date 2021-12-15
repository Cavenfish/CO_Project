from ase.io import read
from ase.io.trajectory import Trajectory
from ase.geometry.analysis import Analysis

systems = Trajectory('moldyn_clus1_geo-opt_10K_40ps_Langevin.traj')

system = systems[0]
system.set_cell([100,100,100])
ana = Analysis(system)
rdf_CC, distances_CC = ana.get_rdf(rmax=8, nbins=160, elements=['C','C'], return_dists=True)[0]
rdf_OO, distances_OO = ana.get_rdf(rmax=8, nbins=160, elements=['O','O'], return_dists=True)[0]

for i in range(1,len(systems)):
	system = systems[i]
	system.set_cell([100,100,100])
	ana = Analysis(system)
	rdf_CC_i, distances_CC_i = ana.get_rdf(rmax=8, nbins=160, elements=['C','C'], return_dists=True)[0]
	rdf_OO_i, distances_OO_i = ana.get_rdf(rmax=8, nbins=160, elements=['O','O'], return_dists=True)[0]
	rdf_CC = [x + y for x, y in zip(rdf_CC, rdf_CC_i)]
	rdf_OO = [x + y for x, y in zip(rdf_OO, rdf_OO_i)]
	print(rdf_CC, distances_CC)



print(rdf_CC, rdf_OO)
rdf_CC_end = []
rdf_OO_end = []
N = len(systems)
for i in range(len(rdf_CC)):
	rdf_CC_end.append(rdf_CC[i]/N)
	rdf_OO_end.append(rdf_OO[i]/N)
	print(rdf_CC_end)
	
with open("Test_rdfs.dat", "a") as f:
	print(rdf_CC_end, rdf_OO_end, file=f)
	print(len(distances_CC), len(rdf_CC_end), file=f)

import matplotlib.pyplot as plt

plt.plot(distances_CC, rdf_CC_end, label='CC')
plt.plot(distances_OO, rdf_OO_end, label='OO')
plt.xlabel('r(Angstrom)')
plt.legend()
plt.show()
