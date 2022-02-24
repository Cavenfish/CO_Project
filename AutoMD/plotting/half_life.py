from ..config import *

def half_life(csvDir, n=51):
    def get_avg():
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

    def plot(avg, ext, property):
        #Half Life formula
        def f(x, E, tau):
            return E * np.exp(-x / tau)
        def g(x, E, r):
            return E - (r*x)

        saveName = csvDir + ext
        x        = avg['Time']/10
        y        = avg[property]
        y2       = savgol_filter(y, n, 2)
        
        p0 = [np.average(y2[0:100]), 1]

        #Run curve fit
        popt, pcov = curve_fit(f, x, y2, p0)
        popt2,pcov2 = curve_fit(g, x, y2, p0)

        #Make plot text
        t1 = popt[1] * np.log(2)
        t2 = (popt2[0]/2) / popt2[1]  
        s  = r'   Exp: $t_{1/2}$ = ' + str(t1)[:5] + ' ps\n' + \
             r'Linear: $t_{1/2}$ = ' + str(t2)[:5] + ' ps'

        plt.plot(x,  y, label='Raw Data')
        plt.plot(x, y2, label='Savitzky-Golay Filter')
        plt.plot(x, f(x, *popt), label='Exponential Fit')
        plt.plot(x, g(x, *popt2), label='Linear Fit')
        plt.text(min(x), min(y), s)
        plt.xlabel('Time (ps)')
        plt.ylabel('Energy (eV)')
        plt.title('Vibrational Energy Dissipation')
        plt.legend()
        plt.tight_layout()
        plt.savefig(saveName)
        plt.close()
        return

        return

    avg = get_avg()
    plot(avg,  'HF_total.png',  'Total Energy')
    plot(avg, 'HF_sliced.png', 'Sliced Energy')

    return
