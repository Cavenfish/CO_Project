from ..config import *

def local_geo(xyz, mol, low, hi):
    def get_rList(atoms):
        pos    = atoms.get_positions()
        masses = atoms.get_masses()
        eCoM   = CoM(pos[mol], masses[mol])
        rList  = []
        for i in range(0, len(masses)//2):
            a     = i*2
            b     = a+2
            cm    = CoM(pos[a:b], masses[a:b])
            diff  = cm - eCoM
            r     = np.sqrt(np.dot(diff,diff))
            rList.append(r)
            rList.append(r)

        return rList

    system = read(xyz)
    rList  = np.array(get_rList(system))
    small  = system[(low <= rList) & (rList <= hi)]

    return small
