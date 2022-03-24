from ..config import *

def diffusion_length(trajFile, saveName, ts=0.1):
    def get_CoM(i,j):
        a   = j*2
        b   = a+2
        pos = traj[i].get_positions()[a:b]
        mas = traj[i].get_masses()[a:b]
        com = CoM(pos, mas)
        return com

    traj = Trajectory(trajFile)
    N    = len(traj)
    n    = len(traj[0]) // 2
    difL = np.zeros(n)

    for i in range(N-1):
        for j in range(n):
            comA     = get_CoM(i  ,j)
            comB     = get_CoM(i+1,j)
            r        = comB - comA
            difL[j] += np.linalg.norm(r)

    t  = N * ts
    df = pd.DataFrame({'Diffusion Length': difL, 'dL/dt': difL/t})
    df.to_csv(saveName)
    return
