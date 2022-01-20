from ase.io import read
from ase.visualize import view

traj = read('pentamer_opt.traj')
view(traj)
