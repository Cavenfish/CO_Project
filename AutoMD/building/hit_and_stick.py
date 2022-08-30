from ..config                import *
from alphashape              import alphashape
from trimesh.proximity       import signed_distance
from scipy.spatial.transform import Rotation

def hit_and_stick(xyz, mol, n, saveName, fmax=1e-6, NVEtime=5, KE=0.25):
    def randVector():
        R     = 1
        theta = np.random.uniform(0, 1) * np.pi
        phi   = np.random.uniform(0, 2) * np.pi

        x = R * np.cos(phi) * np.sin(theta)
        y = R * np.sin(phi) * np.sin(theta)
        z = R * np.cos(theta)

        return [x,y,z]

    def randRotate(v):
        rot = Rotation.random()
        V   = rot.apply(v)
        return V

    def checkDistance(v, pos, r):
        for p in pos:
            a = np.linalg.norm(v[0]-p)
            b = np.linalg.norm(v[1]-p)
            if (a < 3) or (b < 3):
                v[0] += r
                v[1] += r
                v  = checkDistance(v, pos, r)
                break
        return v


    system, _ = prep_system(xyz)
    molec,  _ = prep_system(mol)
    blank_pos = molec.get_positions()
    trajName  = xyz.replace('.xyz', '.traj')

    for i in range(n):
        pos = system.get_positions()
        mas = system.get_masses()
        com = CoM(pos, mas)

        #Get alpha shape
        a = alphashape(pos)

        #Get random unit vector
        r = randVector()
        e = r / np.linalg.norm(r)
        R = 8 * e + com

        #Make new molecule
        new_pos = R + blank_pos
        new_pos = randRotate(new_pos)
        new_pos = checkDistance(new_pos, pos, e)
    
        #Get rho
        diff    = com - CoM(new_pos, mas[:2])
        e       = diff / np.linalg.norm(diff)
        p       = KE * e

        new_rho = [p, p]
        new_mol = Atoms('CO', positions=new_pos, momenta=new_rho)

        #Make system
        system += new_mol

        #Run NVE
        dyn  = VelocityVerlet(system, 1 * units.fs)
        if i == 0:
            traj = Trajectory(trajName, 'w', system)
        else:
            traj = Trajectory(trajName, 'a', system)
        dyn.attach(traj.write, interval=1)
        dyn.run(NVEtime * 1000)

        #Get last frame of simulation
        traj   = Trajectory(trajName)
        system = traj[-1]
        calc   = MvH_CO(atoms=system)
        system.set_calculator(calc)

        #Run geometry optimization
        #traj = Trajectory(trajName, 'a', system)
        opt  = BFGS(system)#, trajectory=traj)
        opt.run(fmax=fmax)

        #Get last frame of simulation
        #system = opt
        #calc   = MvH_CO(atoms=system)
        #system.set_calculator(calc)
        system.set_velocities(np.zeros_like(system.get_velocities()))

    #Write final xyz file
    write(saveName, system)
    return
