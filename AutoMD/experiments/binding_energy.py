from ..config import *
from alphashape import alphashape

def binding_energy(clusters, molecule, n, fmax, saveName,
                   isoClu=None, isoMol=None, alpha=None):
    def randVector():
        R     = 1
        theta = np.random.uniform(0, 1) * np.pi
        phi   = np.random.uniform(0, 2) * np.pi

        x = R * np.cos(phi) * np.sin(theta)
        y = R * np.sin(phi) * np.sin(theta)
        z = R * np.cos(theta)

        return [x,y,z]

    def get_mesh(atoms):
        cmList = []
        for i in range(len(atoms)//2):
            pos = atoms.get_positions()[i*2:i*2+2]
            mas = atoms.get_masses()[i*2:i*2+2]
            com = CoM(pos, mas)
            cmList.append(list(com))

        mesh = alphashape(cmList, alpha)
        return mesh

    def make_name(i, cluster):
        if (isoMol) and (isoClu):
            s    = '_Mol' + str(int(sum(isoMol))) + '_Clu' + \
                   str(int(sum(isoClu)))  + '_BE' + str(i) + '.xyz'
            name = cluster.replace('.xyz', s)
            return name

        if (isoMol):
            s    = '_Mol' + str(int(sum(isoMol))) + '_BE' + str(i) + '.xyz'
            name = cluster.replace('.xyz', s)
            return name

        if (isoClu):
            s    = '_Clu' + str(int(sum(isoClu))) + '_BE' + str(i) + '.xyz'
            name = cluster.replace('.xyz', s)
            return name

        s    = '_BE' + str(i) + '.xyz'
        name = cluster.replace('.xyz', s)

        return name

    #Get all single molecule information
    mol, _    = prep_system(molecule)
    mol_opt   = BFGS(mol)
    mol_opt.run(fmax=fmax)
    mol       = mol_opt.atoms
    mol_E     = mol.get_potential_energy()
    blank_pos = mol.get_positions()
    BE_dict   = {}

    k = 0
    for cluster in clusters:
        k += 1

        #Prep systems
        clu, _    = prep_system(cluster)
        if isoClu:
            m = isoClu * (len(mol.get_masses()) // 2)
            clu.set_masses(m)
        clu_com   = CoM(clu.get_positions(), clu.get_masses())

        #Optimize cluster
        clu_opt   = BFGS(clu)
        clu_opt.run(fmax=fmax)

        #Set optimized structures and get energy
        clu       = clu_opt.atoms
        clu_E     = clu.get_potential_energy()

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

        #Declare a BE dictionary
        cluKey          = 'Cluster ' + str(k)
        BE_dict[cluKey] = []

        for i in range(n):

            #Get random face normal vector
            u = start + i * step
            r = norms[u]
            R = (3 * r) + verts[u]

            #Make new molecule
            new_pos = R + blank_pos
            new_mol = Atoms('CO', positions=new_pos, masses=isoMol)

            #Make system and prep system
            system = clu + new_mol
            calc   = MvH_CO(atoms=system)
            system.set_calculator(calc)

            #Write xyz of binding site pre optimization
            new_name = make_name(i, cluster)
            write(new_name, system)

            #Optimize geometry of new system
            opt = BFGS(system)
            try:
                opt.run(fmax=fmax)
            except:
                BE_dict[cluKey].append(np.nan)
                continue

            #Get energy and binding energy
            ful_E = opt.atoms.get_potential_energy()
            BE    = ful_E - clu_E - mol_E

            #Add energy to dict
            BE_dict[cluKey].append(BE)

    #Make DataFrame and csv
    df = pd.DataFrame(BE_dict)
    df.to_csv(saveName)

    return
