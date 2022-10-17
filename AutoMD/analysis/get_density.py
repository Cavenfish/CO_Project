from ..config import *

def get_density(trajFile, r, center=[]):
    traj = Trajectory(trajFile)
    
    
    v = (4/3) * np.pi * (r**3)
    d = []

    for t in traj:
        if center:
            n = len([x for x in t if np.linalg.norm(x.position -  center) < r])
        else:
            n = len([x for x in t if np.linalg.norm(x.position) < r])
        
        d.append(n/v)

    return d
