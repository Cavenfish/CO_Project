from ..config import *
from ast      import literal_eval

def radial_energy(csvDir, trajDir):
    def plot(df, title, labels, saveName):
        df.plot('t', ['r1', 'r2', 'r3'], label=labels)
        plt.title(title, fontsize=20)
        plt.xlabel('Time (ps)',   fontsize=15)
        plt.ylabel('Energy (eV)', fontsize=15)
        plt.legend(fontsize=10)
        plt.tight_layout()
        plt.savefig(csvDir + saveName)
        plt.close()
        return

    def get_rList(atoms):
        pos    = atoms.get_positions()
        masses = atoms.arrays['masses']
        eCoM   = CoM(pos[:2], masses[:2])
        rList  = []
        for i in range(1, len(masses)//2):
            a     = i*2
            b     = a+2
            cm    = CoM(pos[a:b], masses[a:b])
            diff  = cm - eCoM
            r     = np.sqrt(np.dot(diff,diff))
            rList.append(r)

        return rList
    
    def prep_data(prop):
        def get_r(low, high):
            r = sum(EList[(low < rList) & (rList <  high)])
            try:
                r /= len(rList[(low < rList) & (rList <  high)])
            except:
                r = 0

            return r

        flag = True
        tmp  = {'t':[], 'r1':[], 'r2':[], 'r3':[]}
        N    = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                continue

            try:
                df   = pd.read_csv(csvDir + fname)
                traj = Trajectory(trajDir + fname.replace('.csv', '.traj'))
            except:
                print('Could not find CSV and traj files')
                continue
            
            if len(df) != len(traj):
                print('DF length not equal Traj length')
                continue
            else:
                n = len(df)

            for i in range(n):
                time  = df['Time'][i]
                EList = np.array(literal_eval(df[prop][i]))
                rList = np.array(get_rList(traj[i]))
                
                r1  = get_r(0.0,  4.5)
                r2  = get_r(4.5,  8.0)
                r3  = get_r(8.0, 11.0)

                if flag:
                    tmp['t' ].append(time)
                    tmp['r1'].append(r1)
                    tmp['r2'].append(r2)
                    tmp['r3'].append(r3)
                else:
                    tmp['r1'][i] += r1
                    tmp['r2'][i] += r2
                    tmp['r3'][i] += r3
            N   += 1
            flag = False
        
        tmp['r1'] = np.array(tmp['r1']) / N
        tmp['r2'] = np.array(tmp['r2']) / N
        tmp['r3'] = np.array(tmp['r3']) / N
        
        df = pd.DataFrame(tmp)
        return df
    
    labels = [r'$r_1$', r'$r_2$', r'$r_3$']
    df     = prep_data('All Vibra Energy')
    plot(df, 'Radial Vibrational Energy', labels, 'radVibra.png')

    df     = prep_data('All Trans Energy')
    plot(df, 'Radial Translation Energy', labels, 'radTrans.png')

    df     = prep_data('All Rotat Energy')
    plot(df, 'Radial Rotational Energy', labels, 'radRotat.png')
    return
