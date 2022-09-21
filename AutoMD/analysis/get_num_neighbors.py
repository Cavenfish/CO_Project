from ..config import *

def get_num_neighbors(xyz, mol_index):
    
    def dist(v1,v2):
        d = np.linalg.norm(v1 - v2)
        return d

    atoms = read(xyz)
    coms  = []
    pos   = atoms.get_positions()
    mas   = atoms.get_masses()
    for i in range(len(atoms)//2):
        a = i*2
        b = a+2
        coms.append(CoM(pos[a:b], mas[a:b]))
    
    m     = coms[mol_index]
    tmp   = [x for x in coms if (dist(x,m) <= 5) and (dist(x,m) > 0.5)]
    return len(tmp)
