from ..config import *

def unpickle(file, save, step):
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

    traj  = Trajectory(file)
    # atoms = traj[0]
    # rList = np.array(get_rList(atoms))
    #
    # tmp   = atoms[(low <= rList) & (rList <= hi)]
    # save  = file.replace('.traj', '.xyz')
    # write(save, tmp)

    for frame in traj[1::step]:
        tmp = frame#[(low <= rList) & (rList <= hi)]
        write(save, tmp, append=True)

    return
