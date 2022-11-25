from ..config import *

def all_nu_dissipation(root, noBkg=False):
    def get_avg(csvDir):
        keys = ['Time', 'Sliced Energy']
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
        keys = ['Sliced Energy']
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
        colors   = {'surf':(1,0,1,1), 'subsurf':(0.27,0.04,0.16,1)}
        isotope  = {'co':'CO', '13co':r'$^{13}$CO',
                    '13c18o':r'$^{13}$C$^{18}$O', 'c18o':r'C$^{18}$O'}

        for iso in isotope.keys():
            for dir in os.listdir(root):
                if not os.path.isdir(root + dir):
                    continue

                if dir in colors:
                    labels  = {'surf': 'Surface', 'subsurf': 'Subsurface'}
                    saveName = root + iso + '_' + ext
                    subRoot = root + dir + '/' + iso + '/'
                    for dir2 in os.listdir(subRoot):
                        csvDir  = subRoot + dir2 + '/'
                        if not os.path.isdir(csvDir):
                            continue
                        avg     = get_avg(csvDir)
                        x       = avg['Time']
                        y       = avg[property]
                        l       = labels[dir]
                        c       = colors[dir]
                        _, labs = plt.gca().get_legend_handles_labels()
                        plt.plot(x,  y, color=c,
                                 label=l if l not in labs else "")
                else:
                    continue
                    # csvDir  = root + dir + '/'
                    # avg     = get_avg(csvDir)
                    # x       = avg['Time']
                    # y       = avg[property]
                    # l       = labels[dir]
                    # plt.plot(x,  y, label=l)

            plt.xlabel('Time (ps)',   fontsize=18)
            plt.ylabel('Energy (eV)', fontsize=18)
            plt.title('Vibrational Energy Dissipation', fontsize=20)
            plt.legend(fontsize=12)
            #plt.gca().get_legend().remove()
            plt.tight_layout()
            plt.savefig(saveName, transparent=noBkg)
            plt.close()
        return

    plot('allNu.png', 'Sliced Energy')
    return
