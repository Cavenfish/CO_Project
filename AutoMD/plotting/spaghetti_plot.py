from ..config import *

def spaghetti_plot(csvDir, n=51, eType='Vibrational', noBkg=False):
    def loop(property, ext):
        N     = len(os.listdir(csvDir))
        y_avg = 0
        c     = 255
        k     = 220 / N
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                N -= 1
                continue
            df = pd.read_csv(csvDir + fname)
            x  = df['Time']
            y  = savgol_filter(df[property], n, 2)
            plt.plot(x, y, color=(c/255,0,c/255))
            c -= k

            y_avg += y

        y_avg /= N
        plt.plot(x, y_avg, '--', label='Average', color='black', lw=2)
        plt.ylabel('Energy (eV)', fontsize=18)
        plt.xlabel('Time (ps)',   fontsize=18)
        plt.title(eType + ' Energy Dissipation', fontsize=20)
        plt.legend(fontsize=12)
        plt.savefig(csvDir + ext, transparent=noBkg)
        plt.close()
        return

    loop('Sliced Energy', 'spaghetti_plot.png')
    return
