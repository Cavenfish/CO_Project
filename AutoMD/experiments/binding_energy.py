from ..config                import *
from alphashape              import alphashape
from scipy.spatial.transform import Rotation

def binding_energy(clusters, molecule, n, fmax, saveName, alpha=None):

    def randRotate(v):
        rot = Rotation.random()
        V   = rot.apply(v)
        return V

    def checkDistance(v, pos, r):
        for p in pos:
            a = np.linalg.norm(v[0]-p)
            b = np.linalg.norm(v[1]-p)
            if (a < 3) or (b < 3):
                v[0] += 0.1*r
                v[1] += 0.1*r
                v  = checkDistance(v, pos, r)
                break
        return v
    
    def get_minD(v, pos):
        x = 10000
        for p in pos:
            a = np.linalg.norm(v[0]-p)
            b = np.linalg.norm(v[1]-p)
            if (a < x) or (b < x):
                x = min([a,b])
        return x


    def get_mesh(atoms):
        cmList = []
        for i in range(len(atoms)//2):
            pos = atoms.get_positions()[i*2:i*2+2]
            mas = atoms.get_masses()[i*2:i*2+2]
            com = CoM(pos, mas)
            cmList.append(list(com))

        mesh = alphashape(cmList, alpha)
        return mesh

    #Get all single molecule information
    mol, calc = prep_system(molecule)
    mol_opt   = BFGS(mol)
    mol_opt.run(fmax=fmax)
    mol       = mol_opt.atoms
    mol_E     = mol.get_potential_energy()
    mol_cont  = np.array(calc.get_energy_contributions()[2:])
    blank_pos = mol.get_positions()
    BE_dict   = {}
    contri    = {'Exchange': [], 'Dispersion': [], 'Electrostatic': []}

    k = 0
    for cluster in clusters:
        k += 1

        try:
            clu, calc = prep_system(cluster.replace('.xyz', '_opt.xyz'))
        except:
            clu, calc = prep_system(cluster)
            clu_opt   = BFGS(clu)
            clu_opt.run(fmax=fmax)
            write(cluster.replace('.xyz', '_opt.xyz'), clu)

        #Prep systems
        clu_com   = CoM(clu.get_positions(), clu.get_masses())

        #Set optimized structures and get energy
        clu_E     = clu.get_potential_energy()
        clu_cont  = np.array(calc.get_energy_contributions()[2:])

        #Get variables needed later in loop
        com       = CoM(clu.get_positions(), clu.get_masses())
        mesh      = alphashape(clu.get_positions(), alpha)
        #mesh      = get_mesh(clu)
        norms     = mesh.vertex_normals
        verts     = mesh.vertices
        N         = len(norms)

        #Define method for iterating through possible binding sites
        start     = (N % n) // 2
        step      = N // n

        #Declare a BE and contributions dictionary
        cluKey          = 'Cluster ' + str(k)
        BE_dict[cluKey] = []

        for i in range(n):

            #Check if calculation already done
            s     = '_BE%d_opt.xyz' % i
            check = cluster.replace('.xyz', s)
            try:
                sys, calc = prep_system(check)
                ful_E     = sys.get_potential_energy()
                ful_cont  = np.array(calc.get_energy_contributions()[2:])
                BE        = ful_E - clu_E - mol_E
                cont      = ful_cont - clu_cont - mol_cont
                cont      = np.abs(cont) / sum(np.abs(cont)) * 100
                BE_dict[cluKey].append(BE)
                contri['Exchange'].append(cont[0])
                contri['Dispersion'].append(cont[1])
                contri['Electrostatic'].append(cont[2])
                continue
            except:
                pass

            #Get random face normal vector
            u = start + i * step
            r = norms[u]
            R = 0.1*r + verts[u]

            #Make new molecule
            tmp_pos = randRotate(blank_pos)
            new_pos = R + tmp_pos
            new_pos = checkDistance(new_pos, clu.get_positions(), r)
            new_mol = Atoms('CO', positions=new_pos)

            D = get_minD(new_pos, clu.get_positions())
            if D > 5:
                BE_dict[cluKey].append(np.nan)
                contri['Exchange'].append(np.nan)
                contri['Dispersion'].append(np.nan)
                contri['Electrostatic'].append(np.nan)
                continue

            #Make system and prep system
            system = clu + new_mol
            calc   = MvH_CO(atoms=system)
            system.set_calculator(calc)

            #Write xyz of binding site pre optimization
            s        = '_BE' + str(i) + '.xyz'
            new_name = cluster.replace('.xyz', s)
            write(new_name, system)

            #Optimize geometry of new system
            opt = BFGS(system)
            try:
                opt.run(fmax=fmax, steps=5000)
                assert opt.converged() == True
            except:
                BE_dict[cluKey].append(np.nan)
                contri['Exchange'].append(np.nan)
                contri['Dispersion'].append(np.nan)
                contri['Electrostatic'].append(np.nan)
                continue

            #Write out opt geo
            opt_name = new_name.replace('.xyz', '_opt.xyz')
            write(opt_name, system)

            #Get energy and binding energy
            ful_E = opt.atoms.get_potential_energy()
            BE    = ful_E - clu_E - mol_E
    
            #Get energy contributions
            ful_cont  = np.array(calc.get_energy_contributions()[2:])
            cont      = ful_cont - clu_cont - mol_cont
            cont      = np.abs(cont) / sum(np.abs(cont)) * 100
            
            #Add energy contributions to dicts
            contri['Exchange'].append(cont[0])
            contri['Dispersion'].append(cont[1])
            contri['Electrostatic'].append(cont[2])

            #Add energy to dict
            BE_dict[cluKey].append(BE)


    #Make DataFrame and csv
    df = pd.DataFrame(BE_dict)
    df.to_csv(saveName)
    df = pd.DataFrame(contri)
    df.to_csv(saveName.replace('.csv', '_contri.csv'))

    return
