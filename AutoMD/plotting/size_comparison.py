from ..config import *
from operator import itemgetter

def size_comparison(root, noBkg=False):
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

    def get_std(csvDir, avg):
        keys = ['Total Energy', 'Sliced Energy']
        tmp  = {}
        N    = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                continue
            N += 1
            df = pd.read_csv(csvDir + fname)

            for key in keys:
                if key not in tmp:
                    tmp[key]  = (df[key] - avg[key])**2
                else:
                    tmp[key] += (df[key] - avg[key])**2
        for key in tmp:
            tmp[key] /= (N - 1)
            tmp[key]  = np.sqrt(tmp[key])

        std = pd.DataFrame(tmp)
        return std

    def plot(ext, property):
        saveName = root + ext
        labels   = {'60co': '60 CO', '160co':'160 CO', '200co': '200 CO',
                    '240co': '240 CO'}
        colors   = {'60co': (1,0,1,1), '160co':(0,0.75,0.75,1),
                    '200co': (0.27,0.04,0.16,1), '240co': (0.94,0.31,0.15,1)}
        for dir in os.listdir(root):
            if not os.path.isdir(root + dir):
                continue

            csvDir  = root + dir + '/'
            avg     = get_avg(csvDir)
            #std     = get_std(csvDir, avg)
            x       = avg['Time']
            y       = avg[property]
            #yerr    = std[property]
            l       = labels[dir]
            c       = colors[dir]
            plt.plot(x,  y, label=l, color=c)
            #mks, cs, bs = plt.errorbar(x,  y, yerr, label=l, color=c,
            #                           elinewidth=0.5, capsize=0.75)

            #[bar.set_alpha(0.05) for bar in bs]
            #[cap.set_alpha(0.05) for cap in cs]

        plt.xlabel('Time (ps)',   fontsize=18)
        plt.ylabel('Energy (eV)', fontsize=18)
        plt.title('Vibrational Energy Dissipation', fontsize=20)

        h, l = plt.gca().get_legend_handles_labels()
        tmp  = sorted(zip(h,l), key=itemgetter(1))
        h, l = zip(*tmp)

        plt.legend(h, l, fontsize=12)
        plt.tight_layout()
        plt.savefig(saveName, transparent=noBkg)
        plt.close()
        return

    plot('sizeComparison.png', 'Total Energy')
    plot('sizeComparison_sliced.png', 'Sliced Energy')
    return
