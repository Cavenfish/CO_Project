import gc, sys
from ..config                import *
from scipy.spatial.transform import Rotation
from yaml                    import safe_load

def hit_and_stick(inp):
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
            if (a < 5) or (b < 5):
                v[0] += 0.1*r
                v[1] += 0.1*r
                v  = checkDistance(v, pos, r)
                break
        return v

    def convertUnits(KE):
        E = float(KE.split('_')[0])
        u = KE.split('_')[1]

        if u == 'K':
            return E * 8.6173e-5
        elif u == 'kJ':
            return E * 6.2415e21
        elif u == 'eV':
            return E
        else:
            sys.exit(f'Wrong Energy Units\n{u} is not allowed')

    #Read in input card
    with open(inp, 'r') as f:
        p = safe_load(f)

    #Initialize system 
    n         = p['size'] - 1
    system    = read(p['xyz'])
    molec     = read(p['mol'])
    blank_pos = molec.get_positions()
    coMass    = molec.get_masses()

    #Make blank traj file 
    if 'trajName' in p:
        interval = p['trajInter']
        with open(p['trajName'], 'w') as f:
            pass

    #Run hit and stick
    for i in range(n):
        pos = system.get_positions()
        mas = system.get_masses()
        com = CoM(pos, mas)

        #Get random unit vector
        r = randVector()
        e = r / np.linalg.norm(r)
        R = 8 * e + com

        #Make new molecule
        tmp_pos = randRotate(blank_pos)
        new_pos = R + tmp_pos
        new_pos = checkDistance(new_pos, pos, e)
    
        #Get rho unit vector
        diff    = com - CoM(new_pos, mas[:2])
        e       = diff / np.linalg.norm(diff)
        
        #Get KE in eV
        KE = convertUnits(p['NVE']['KE'])
        
        #Make rho
        mv   = np.sqrt(KE * 2 * coMass)
        rho0 = mv[0] * e
        rho1 = mv[1] * e

        #Prep incoming molecule 
        new_rho = [rho0, rho1]
        new_mol = Atoms('CO', positions=new_pos, momenta=new_rho)

        #Make system
        system += new_mol
        calc    = MvH_CO(atoms=system)
        system.set_calculator(calc)

        #Prep trajectory file
        if 'trajName' in p:
            traj = Trajectory(p['trajName'], 'a', system)

        #Run NVE
        if 'NVE' in p:
            dyn = VelocityVerlet(system, 1 * units.fs)
            if 'trajName' in p: dyn.attach(traj.write, interval=interval)
            dyn.run(p['NVE']['time'] * 1000)

        #Run Geo Opt
        if 'OPT' in p:
            fmax = float(p['OPT']['fmax']) 
            opt  = BFGS(system)
            if 'trajName' in p: opt.attach(traj.write, interval=interval)
            opt.run(fmax=fmax)

        #Run NVT
        if 'NVT' in p:
            T   = p['NVT']['T']
            mu  = p['NVT']['mu']
            dyn = Langevin(system, 1 * units.fs, friction=mu, temperature_K=T)
            if 'trajName' in p: dyn.attach(traj.write, interval=interval)
            dyn.run(p['NVT']['time'] * 1000)

        #Collect garbage (might help with performace)
        gc.collect()

    #Write final xyz file
    write(p['saveName'], system)
    return
