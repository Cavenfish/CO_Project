from ..config import *
from alphashape import alphashape
import warnings

def binding_energy(clusters, molecule, n, fmax, alpha=None):
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

        #Prevent overasking for binding sites
        if n > N:
            warnings.warn('\n Requested number of binding sites is larger than ' + \
                          'number of verticies. The requested amount will ' + \
                          'be made equal to the number of verticies.')
            continue



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
            new_mol = Atoms('CO', positions=new_pos)

            #Make system and prep system
            system = clu + new_mol
            calc   = MvH_CO(atoms=system)
            system.set_calculator(calc)

            #Write xyz of binding site pre optimization
            new_name = cluster.replace('.xyz', '_BE' + str(i) + '.xyz')
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
    df.to_csv('BE_test.csv')

    return
