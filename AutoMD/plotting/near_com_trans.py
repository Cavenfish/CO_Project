from ..config import *

def near_com_trans(csvDir, low, hi):
    def get_rList(atoms):
        pos    = atoms.get_positions()
        masses = atoms.get_masses()
        eCoM   = CoM(pos[:2], masses[:2])
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
    return
