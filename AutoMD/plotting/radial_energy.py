from ..config import *
from ast      import literal_eval

def radial_energy(csvDir, trajDir):
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
                time  = i * 0.1 # in ps 
                EList = np.array(literal_eval(df[prop][i]))
                rList = np.array(get_rList(traj[i]))
                
                r1  = sum(EList[(2.5 < rList) & (rList <  4.5)])
                r2  = sum(EList[(6.0 < rList) & (rList <  8.0)])
                r3  = sum(EList[(9.0 < rList) & (rList < 11.0)])

                r1 /= len(rList[(2.5 < rList) & (rList <  4.5)])
                r2 /= len(rList[(6.0 < rList) & (rList <  8.0)])
                r3 /= len(rList[(9.0 < rList) & (rList < 11.0)])

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
    
    df = prep_data('All Vibra Energy')
    df.plot('t', ['r1', 'r2', 'r3'])
    plt.show()
    return
