from ..config import *
from yaml import dump

def get_radial_energy(traj, saveName=False):
    def Dist(c1,c2):
        d = np.linalg.norm(c2-c1)
        return d

    saveName = saveName if saveName else traj.replace('.traj', '.yaml')

    tj = Trajectory(traj)
    l  = {}

    for i in range(len(tj)):
        x   = tj[i]
        num = len(x)
        tim = i * 0.1
        pos = x.positions
        mas = x.get_masses()
        vel = x.get_velocities()
        com = CoM(pos[0:2], mas[0:2])
        tmp = {'dist':[], 'E_vib':[], 'E_rot':[], 'E_tra':[]}

        #Get lists
        for j in range(0,num,2):
            tmp['dist'].append(Dist(com, CoM(pos[j:j+2], mas[j:j+2])))
            tmp['E_vib'].append(calc_Evib(pos[j:j+2],mas[j:j+2],vel[j:j+2]))
            tmp['E_rot'].append(calc_Erot(pos[j:j+2],mas[j:j+2],vel[j:j+2]))
            tmp['E_tra'].append(calc_Etran(pos[j:j+2],mas[j:j+2],vel[j:j+2]))

        l[tim] = tmp

    with open(saveName, 'w') as f:
        dump(l, f)

    return l
