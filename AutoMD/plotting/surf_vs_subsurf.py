from ..config import *

def surf_vs_subsurf(root, noBkg=False):
    def get_avg(csvDir):
        keys = ['Time', 'Total Energy', 'Sliced Energy']
        tmp  = {}
        N    = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                continue
            N += 1
            df = pd.read_csv(csvDir + fname)

            for key in keys:
                if key not in tmp:
                    tmp[key]  = df[key]
                else:
                    tmp[key] += df[key]
        for key in tmp:
            tmp[key] /= N

        avg = pd.DataFrame(tmp)
        return avg

    def get_dist_range(csvDir):
         
        return

    def plot(ext, property):
        saveName = root + ext
        labels   = {'surf':'Surface', 'subsurf':'Sub-Surface'}
        colors   = {'surf':(1,0,1,1), 'subsurf':(0.27,0.04,0.16,1)}
        for dir in os.listdir(root):
            if not os.path.isdir(root + dir):
                continue

            csvDir  = root + dir + '/'
            avg     = get_avg(csvDir)
            x       = avg['Time']
            y       = avg[property]
            l       = labels[dir]
            c       = colors[dir]
            plt.plot(x,  y, label=l, color=c)

        plt.xlabel('Time (ps)',   fontsize=18)
        plt.ylabel('Energy (eV)', fontsize=18)
        plt.title('Vibrational Energy Dissipation', fontsize=20)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(saveName, transparent=noBkg)
        plt.close()
        return

    plot('surfVsSub.png', 'Total Energy')
    plot('surfVsSub_sliced.png', 'Sliced Energy')
    return
