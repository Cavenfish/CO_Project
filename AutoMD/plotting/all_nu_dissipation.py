from ..config import *

def all_nu_dissipation(root, noBkg=False):
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
        labels   = {'1': r'$\nu$ = 1', '2': r'$\nu$ = 2', '3': r'$\nu$ = 3',
                    '4': r'$\nu$ = 4', '5': r'$\nu$ = 5', '6': r'$\nu$ = 6',
                'nu1': r'$\nu$ = 1', 'nu2': r'$\nu$ = 2', 'nu3': r'$\nu$ = 3',
                'nu4': r'$\nu$ = 4', 'nu5': r'$\nu$ = 5', 'nu6': r'$\nu$ = 6'}
        #colors   = {'surf':(1,0,1,1), 'subsurf':(0.27,0.04,0.16,1)}
        for dir in os.listdir(root):
            if not os.path.isdir(root + dir):
                continue

            csvDir  = root + dir + '/'
            avg     = get_avg(csvDir)
            x       = avg['Time']
            y       = avg[property]
            l       = labels[dir]
            #c       = colors[dir]
            plt.plot(x,  y, label=l)#, color=c)

        plt.xlabel('Time (ps)',   fontsize=18)
        plt.ylabel('Energy (eV)', fontsize=18)
        plt.title('Vibrational Energy Dissipation', fontsize=20)
        #plt.legend(fontsize=12)
        #plt.gca().get_legend().remove()
        plt.tight_layout()
        plt.savefig(saveName, transparent=noBkg)
        plt.close()
        return

    plot('allNu.png', 'Total Energy')
    plot('allNu_sliced.png', 'Sliced Energy')
    return
