from ..config import *

def get_num_neighbors(xyz, mol, D):
    
    def dist(v1,v2,v3):
        a = np.linalg.norm(v1 - v2)
        b = np.linalg.norm(v1 - v3)
        d = min([a,b])
        return d

    def find_pair(n):
        if n%2:
            p = (n-1, n)
        else:
            p = (n, n+1)
        return p
    
    atoms = read(xyz)
    tmp   = [x % len(atoms) for x in mol]
    pos   = atoms.get_positions()
    a1    = pos[tmp[0]]
    a2    = pos[tmp[1]]
    pos   = [list(pos[i]) for i in range(len(pos)) if i not in mol]
    tmp   = [pos.index(x) for x in pos if dist(x,a1,a2) <= D]
    tmp   = [find_pair(x) for x in tmp]
    N     = len(set(tmp))
    
    return N
