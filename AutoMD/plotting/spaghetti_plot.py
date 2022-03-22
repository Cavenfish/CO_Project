from ..config import *

def spaghetti_plot(csvDir, n=51):
    def loop(property, ext):
        N     = len(os.listdir(csvDir))
        y_avg = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                N -= 1
                continue
            df = pd.read_csv(csvDir + fname)
            x  = df['Time']
            y  = savgol_filter(df[property], n, 2)
            plt.plot(x, y)

            y_avg += y

        y_avg /= N
        plt.plot(x, y_avg, '--', label='Average', color='black', lw=2)
        plt.ylabel('Energy (eV)')
        plt.xlabel('Time (ps)')
        plt.title('Vibrational Energy Dissipation')
        plt.legend()
        plt.savefig(csvDir + ext)
        plt.close()
        return

    loop( 'Total Energy', 'spaghetti_plot.png')
    loop('Sliced Energy', 'spaghetti_plot_sliced.png')
    return
