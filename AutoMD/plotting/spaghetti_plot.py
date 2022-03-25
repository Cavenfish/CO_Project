from ..config import *

def spaghetti_plot(csvDir, n=51):
    def loop(property, ext):
        N     = len(os.listdir(csvDir))
        y_avg = 0
        c     = 255
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                N -= 1
                continue
            df = pd.read_csv(csvDir + fname)
            x  = df['Time']
            y  = savgol_filter(df[property], n, 2)
            plt.plot(x, y, color=(c/255,0,c/255))
            c -= 1.5

            y_avg += y

        y_avg /= N
        plt.plot(x, y_avg, '--', label='Average', color='black', lw=2)
        plt.ylabel('Energy (eV)', fontsize=15)
        plt.xlabel('Time (ps)',   fontsize=15)
        plt.title('Vibrational Energy Dissipation', fontsize=20)
        plt.legend(fontsize=10)
        plt.savefig(csvDir + ext, transparent=True)
        plt.close()
        return

    loop( 'Total Energy', 'spaghetti_plot.png')
    loop('Sliced Energy', 'spaghetti_plot_sliced.png')
    return
