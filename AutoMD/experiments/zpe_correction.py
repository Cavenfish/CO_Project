from ..config import *

def zpe_correction(molecule, s0, s1, n, m, saveName):

    def get_EsZPE(xyz):
        system, _ = prep_system(xyz)
        E         = system.get_potential_energy()
        vib       = Vibrations(system, delta=0.0001)
        vib.run()
        vib.summary(log=xyz.replace('.xyz', '_vib.out'))
        zpe       = vib.get_zero_point_energy()
        vib.clean()
        EsZPE     = E + zpe
        return EsZPE

    dic = {}
    mol = get_EsZPE(molecule)

    for i in range(n):
        key      = 'Cluster %d' % i
        dic[key] = []
        for j in range(m):
            try:
                clu = get_EsZPE(s0 % i)
                ful = get_EsZPE(s1 % (i,j))
                dic[key].append(ful - clu - mol)
            except:
                dic[key].append(np.nan)

    df  = pd.DataFrame(dic)
    df.to_csv(saveName)
    return
