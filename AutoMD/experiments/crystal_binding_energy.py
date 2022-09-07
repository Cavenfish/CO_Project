from ..config import *
from scipy.spatial.transform import Rotation
from ase.constraints import FixAtoms

def crystal_binding_energy(cluster, molecule, n, posList, fmax, saveName):

    def randRotate(v):
        rot = Rotation.random()
        V   = rot.apply(v)
        return V

    #Get all single molecule information
    mol, _    = prep_system(molecule)
    mol_opt   = BFGS(mol)
    mol_opt.run(fmax=fmax)
    mol       = mol_opt.atoms
    mol_E     = mol.get_potential_energy()
    blank_pos = mol.get_positions()
    cluKey    = 'Binding Energy'
    BE_dict   = {cluKey:[]}

    #Prep systems
    clu, _    = prep_system(cluster)
    fixed     = FixAtoms(
                    indices=[x.index for x in clu if x.position[2] < 32])
    clu.set_constraint(fixed)

    #Optimize cluster
    clu_opt   = BFGS(clu, trajectory=cluster.replace('xyz','traj'))
    clu_opt.run(fmax=fmax)

    #Set optimized structures and get energy
    clu       = clu_opt.atoms
    clu_E     = clu.get_potential_energy()
    (x,y,z)   = clu.get_positions()[-1]

    for i in range(len(posList)):
        for j in range(n):

            #Get new molecule pos
            R       = np.array(posList[i])

            #Make new molecule
            tmp_pos = randRotate(blank_pos)
            new_pos = R + tmp_pos
            new_mol = Atoms('CO', positions=new_pos)

            #Make system and prep system
            system = clu + new_mol
            fixed  = FixAtoms(
                    indices=[x.index for x in system if x.position[2] < 32])
            system.set_constraint(fixed)
            calc   = MvH_CO(atoms=system)
            system.set_calculator(calc)

            #Write xyz of binding site pre optimization
            s        = '_pos%d_BE%d.xyz' % (i,j)
            new_name = cluster.replace('.xyz', s)
            write(new_name, system)

            #Optimize geometry of new system
            opt = BFGS(system, trajectory=s.replace('.xyz', '.traj'))
            try:
                opt.run(fmax=fmax)
            except:
                BE_dict[cluKey].append(np.nan)
                continue

            #Write out opt geo
            opt_name = new_name.replace('.xyz', '_opt.xyz')
            write(opt_name, system)

            #Get energy and binding energy
            ful_E = opt.atoms.get_potential_energy()
            BE    = ful_E - clu_E - mol_E

            #Add energy to dict
            BE_dict[cluKey].append(BE)

    #Make DataFrame and csv
    df = pd.DataFrame(BE_dict)
    df.to_csv(saveName)

    return
