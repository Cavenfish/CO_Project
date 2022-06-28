from ..config import *
import matplotlib as mpl
from ast      import literal_eval

def neighbor_energy(csvDir, trajDir, noBkg=False):
    def plot(df, title, labels, saveName):
        #mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
        #    color=[(1,0,1,1), (0.27,0.04,0.16,1), (0,0.75,0.75,1)])
        df.plot('t', df.keys()[1:])
        plt.title(title, fontsize=20)
        plt.xlabel('Time (ps)',   fontsize=18)
        plt.ylabel('Energy (eV)', fontsize=18)
        #plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(csvDir + saveName, transparent=noBkg)
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

    def plot_data(prop, title, labels, saveName):
        def get_r(low, high):
            r = EList[(low < rList) & (rList <  high)]
            return r

        flag = True
        tmp  = {'t':[]}
        N    = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                continue
            N+=1
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

                tmp['t' ].append(time)
                k = 0
                for j in r1:
                    k  += 1
                    key = 'mol' + str(k)
                    if flag:
                        tmp[key] = [j]
                    else:
                        try:
                            tmp[key].append(j)
                        except:
                            try:
                                tmp[key].append(0)
                            except:
                                k -= 1
                                continue
                flag = False

            del_this = []
            for i in tmp:
                if len(tmp[i]) != len(tmp['t']):
                    del_this.append(i)
            for j in del_this:
                del(tmp[j])
            df = pd.DataFrame(tmp)
            plot(df, title, labels, saveName + str(N) + '.png')
        return df

    labels = [r'$r_1$']
    df     = plot_data('All Vibra Energy', 'Radial Vibrational Energy',
                        labels, 'radVibra')

    df     = prep_data('All Trans Energy', 'Radial Translation Energy',
                       labels, 'radTrans.png')

    df     = prep_data('All Rotat Energy', 'Radial Rotational Energy',
                       labels, 'radRotat.png')
    return
