from ..config import *
from operator import itemgetter

def all_isotopes_dissipation(root):
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
        labels   = {'co':'CO', '13co':r'$^{13}$CO', 
                    '13c18o':r'$^{13}$C$^{18}$O', 'c18o':r'C$^{18}$O'}
        shapes   = {'co':'-', '13co':':', 'c18o':'--', '13c18o':'-.'}
        colors   = {'co':(1,0,1,1), '13co':(0,0.9,0.9,1), 
                    '13c18o':(0.27,0.04,0.16,1), 'c18o':(0.94,0.31,0.15,1)}   
        for dir in os.listdir(root):
            if not os.path.isdir(root + dir):
                continue

            csvDir  = root + dir + '/' 
            avg     = get_avg(csvDir)
            x       = avg['Time']
            y       = avg[property]
            l       = labels[dir]
            c       = colors[dir]
            s       = shapes[dir]
            plt.plot(x,  y, label=l, color=c, ls=s)
        
        plt.xlabel('Time (ps)',   fontsize=15)
        plt.ylabel('Energy (eV)', fontsize=15)
        plt.title('Vibrational Energy Dissipation', fontsize=20)

        h, l = plt.gca().get_legend_handles_labels()
        tmp  = sorted(zip(h,l), key=itemgetter(1))
        h, l = zip(*tmp)
        plt.legend(h, l, fontsize=10)

        plt.tight_layout()
        plt.savefig(saveName)
        plt.close()
        return
    
    plot('all_dissipation.png', 'Total Energy')
    plot('all_dissipation_sliced.png', 'Sliced Energy')
    return
