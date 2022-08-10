from ..config import *

def crystal_binding_energy(cluster, molecule, n, fmax, saveName):

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
    clu_com   = CoM(clu.get_positions(), clu.get_masses())

    #Optimize cluster
    clu_opt   = BFGS(clu)
    clu_opt.run(fmax=fmax)

    #Set optimized structures and get energy
    clu       = clu_opt.atoms
    clu_E     = clu.get_potential_energy()
    (x,y,z)   = clu.get_positions()[-1]

    for i in range(n):
        for j in range(n):

            #Get new molecule pos
            R       = [(i*x/n)+3, (j*y/n), 0]

            #Make new molecule
            new_pos = R + blank_pos
            new_mol = Atoms('CO', positions=new_pos)

            #Make system and prep system
            system = clu + new_mol
            calc   = MvH_CO(atoms=system)
            system.set_calculator(calc)

            #Write xyz of binding site pre optimization
            s        = '_BE%d-%d.xyz' % (i,j)
            new_name = cluster.replace('.xyz', s)
            write(new_name, system)

            #Optimize geometry of new system
            opt = BFGS(system)
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
