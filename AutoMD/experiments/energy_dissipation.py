from ..config import *
from alphashape import alphashape
import random

def energy_dissipation(xyz, loc, nu, n, iso=[], N=500000):
    def label_mols(atoms):
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

        pos  = atoms.get_positions()
        rlp  = range(len(pos))
        alp  = alphashape(pos, None)
        tmp  = [find_pair(i) for i in rlp if np.array(pos[i]) in alp.vertices]
        surf = set(tmp)
        sub  = []
        for p in surf:
            a1    = pos[p[0]]
            a2    = pos[p[1]]
            for i in rlp:
                x = pos[i]
                if dist(x,a1,a2) <= 3.6:
                    sub.append(i)
        sub  = set([find_pair(i) for i in sub])
        tmp  = set([find_pair(i) for i in rlp])
        bulk = [p for p in tmp if (p not in sub) and (p not in surf)]

        dic  = {'surf':surf, 'subsurf':sub, 'bulk':bulk}
        return dic

    atoms  = read(xyz)
    labels = label_mols(atoms)

    swap   = list(random.choice(list(labels[loc])))
    if iso:
        xyz  = add_isotope(xyz, swap, iso)
        swap = [0,1]

    traj   = run_verletMD(xyz, n=20000, i=100)
    xyz    = traj.replace('.traj', '.xyz')
    traj   = Trajectory(traj)
    write(xyz, traj[-1])
    E      = Morse_excitation(nu, n)
    xyz    = excite_molecule(xyz, swap, E)
    traj   = run_verletMD(xyz, n=N, i=100)
    df     = make_NVE_output(traj, traj.replace('.traj', '.csv'), 0.1)
    return
