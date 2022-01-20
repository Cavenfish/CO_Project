import sys, time
from ase.io import read
from ase import units, Atoms
from MvH_CO_JM8 import MvH_CO
from ase.visualize import view
from ase.vibrations import Vibrations
from ase.md.verlet import VelocityVerlet

def view_xyz(file):
    obj = read(file)
    view(file)
    return 


