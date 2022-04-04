from ..config import *
import matplotlib as mpl

def energy_contributions(csvDir, n=51, noBkg=False):
    def get_avg():
        tmp = {}
        N   = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                continue
            N += 1
            df = pd.read_csv(csvDir + fname)

            for key in df:
                if 'All' in key:
                    continue
                if key not in tmp:
                    tmp[key]  = df[key]
                else:
                    tmp[key] += df[key]
        for key in tmp:
            tmp[key] /= N

        avg = pd.DataFrame(tmp)
        return avg

    def plot(avg, ext, label, property, title, all=[]):
        mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
            color=[(1,0,1,1), (0.27,0.04,0.16,1), (0,0.75,0.75,1)])
        saveName = csvDir + ext
        x        = avg['Time']
        y        = savgol_filter(avg[property], n, 2)
        if all:
            y2   = savgol_filter(avg[all[0]], n, 2)
            y3   = savgol_filter(avg[all[2]], n, 2)
            plt.plot(x, y2, label=all[1])
            plt.plot(x, y3, label=all[3])
        plt.plot(x,  y, label=label)
        plt.xlabel('Time (ps)',   fontsize=18)
        plt.ylabel('Energy (eV)', fontsize=18)
        plt.title(title, fontsize=18)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(saveName, transparent=noBkg)
        plt.close()
        return

    avg = get_avg()
    plot(avg, 'avgTrans.png', 'Translational', 'Avg Trans Energy',
         'Average Energy of non-Excited Molecules')

    plot(avg, 'avgRot.png', 'Rotational', 'Avg Rotat Energy',
         'Average Energy of non-Excited Molecules')

    plot(avg, 'avgVibra.png', 'Vibrational', 'Avg Vibra Energy',
         'Average Energy of non-Excited Molecules')

    plot(avg, 'oneTrans.png', 'Translational', 'One Trans Energy',
         'Energy of Excited Molecule')

    plot(avg, 'oneRot.png', 'Rotational', 'One Rotat Energy',
         'Energy of Excited Molecule')

    plot(avg, 'oneVibra.png', 'Vibrational', 'One Vibra Energy',
         'Energy of Excited Molecule')

    plot(avg, 'avg.png', 'Translational', 'Avg Trans Energy',
         'Average Energy of non-Excited Molecules',
         ['Avg Rotat Energy', 'Rotational',
          'Avg Vibra Energy', 'Vibrational'])

    plot(avg, 'one.png', 'Translational', 'One Trans Energy',
         'Energy of Excited Molecule',
         ['One Rotat Energy', 'Rotational',
          'One Vibra Energy', 'Vibrational'])

    return
