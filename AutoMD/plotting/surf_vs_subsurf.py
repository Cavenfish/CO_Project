from ..config import *

def surf_vs_subsurf(root, noBkg=False, errBar=False):
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
        labels   = {'surf':'Surface', 'subsurf':'Sub-Surface'}
        colors   = {'surf':(1,0,1,1), 'subsurf':(0.27,0.04,0.16,1)}
        isotope  = {'co':'CO', '13co':r'$^{13}$CO',
                    '13c18o':r'$^{13}$C$^{18}$O', 'c18o':r'C$^{18}$O'}
        for iso in isotope.keys():
            for loc in labels.keys():
                csvDir  = root + loc + '/' + iso + '/nu1/'
                if not os.path.isdir(csvDir):
                    continue

                avg     = get_avg(csvDir)
                std     = get_std(csvDir, avg)
                x       = avg['Time']
                y       = avg[property]
                l       = labels[loc]
                c       = colors[loc]
                if errBar:
                    yerr    = std[property]
                    mks, cs, bs = plt.errorbar(x,  y, yerr, label=l, color=c,
                                           elinewidth=0.5, capsize=0.75)
                    [bar.set_alpha(0.05) for bar in bs]
                    [cap.set_alpha(0.05) for cap in cs]
                else:
                    plt.plot(x, y, label=l, color=c)

            title    = isotope[iso] + ' Vibrational Energy Dissipation'
            saveName = root + iso + '_' + ext

            plt.xlabel('Time (ps)',   fontsize=18)
            plt.ylabel('Energy (eV)', fontsize=18)
            plt.title(title, fontsize=20)
            plt.legend(fontsize=12)
            plt.tight_layout()
            plt.savefig(saveName, transparent=noBkg)
            plt.close()
        return

    plot('surfVsSub.png', 'Total Energy')
    plot('surfVsSub_sliced.png', 'Sliced Energy')
    return
