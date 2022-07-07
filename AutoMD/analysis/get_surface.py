from ..config import *
from alphashape import alphashape

def get_surface(xyz, alpha):
    atoms  = read(xyz)
    cmList = []
    for i in range(len(atoms)//2):
        pos = atoms.get_positions()[i*2:i*2+2]
        mas = atoms.get_masses()[i*2:i*2+2]
        com = CoM(pos, mas)
        cmList.append(list(com))

    alpha_shape = alphashape(cmList, alpha)
    tmp = [cmList.index(x) for x in cmList if np.array(x) in alpha_shape.vertices]

    surface = Atoms()
    for i in range(len(atoms)):
        if i in tmp:
            surface += atoms[i*2:i*2+2]
        else:
            continue

    return surface
