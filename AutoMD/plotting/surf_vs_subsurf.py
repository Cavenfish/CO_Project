from ..config import *

def surf_vs_subsurf(root):
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

    def plot(ext, property):
        saveName = root + ext
        labels   = {'surf':'Surface', 'subsurf':'Sub-Surface'}
        for dir in os.listdir(root):
            if not os.path.isdir(root + dir):
                continue

            csvDir  = root + dir + '/'
            avg     = get_avg(csvDir)
            x       = avg['Time']
            y       = avg[property]
            plt.plot(x,  y, label=labels[dir])

        plt.xlabel('Time (ps)')
        plt.ylabel('Energy (eV)')
        plt.title('Vibrational Energy Dissipation')
        plt.legend()
        plt.tight_layout()
        plt.savefig(saveName)
        plt.close()
        return

    plot('surfVsSub.png', 'Total Energy')
    plot('surfVsSub_sliced.png', 'Sliced Energy')
    return
