from ..config import *

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
        for dir in os.listdir(root):
            if not os.path.isdir(root + dir):
                continue

            csvDir  = root + dir + '/' 
            avg     = get_avg(csvDir)
            x       = avg['Time']
            y       = avg[property]
            plt.plot(x,  y, label=labels[dir])
        
        plt.xlabel('Time (ps)',   fontsize=15)
        plt.ylabel('Energy (eV)', fontsize=15)
        plt.title('Vibrational Energy Dissipation', fontsize=20)
        plt.legend(fontsize=10)
        plt.tight_layout()
        plt.savefig(saveName)
        plt.close()
        return
    
    plot('all_dissipation.png', 'Total Energy')
    plot('all_dissipation_sliced.png', 'Sliced Energy')
    return
