from ..config          import *
from alphashape        import alphashape
from trimesh.proximity import signed_distance

def hit_and_stick(xyz, n, saveName, fmax=1e-6):
    def randVector():
        R     = 1
        theta = np.random.uniform(0, 1) * np.pi
        phi   = np.random.uniform(0, 2) * np.pi

        x = R * np.cos(phi) * np.sin(theta)
        y = R * np.sin(phi) * np.sin(theta)
        z = R * np.cos(theta)

        return [x,y,z]

    system, _ = prep_system(xyz)
    blank_pos = system.get_positions()
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
        p = 0.25 * -e

        #Check is random spawn spot is far enough
        try:
            d = signed_distance(a, [R])[0]
            x =  6 + d
            R = (8 + x) * e + com
        except:
            pass

        #Make new molecule
        new_pos = R + blank_pos
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
        dyn.run(5000)

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
        traj   = Trajectory(trajName)
        system = traj[-1]
        calc   = MvH_CO(atoms=system)
        system.set_calculator(calc)
        system.set_velocities(np.zeros_like(system.get_velocities()))

    #Write final xyz file
    write(saveName, system)
    return
