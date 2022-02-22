from ..config import *

def spaghetti_plot(csvDir, n=51):

    N     = len(os.listdir(csvDir))
    y_avg = 0
    for fname in os.listdir(csvDir):
        if '.csv' not in fname:
            N -= 1
            continue
        df = pd.read_csv(csvDir + fname)
        x  = df['Time']/10
        y  = savgol_filter(df['Total Energy'], n, 2)
        plt.plot(x, y)

        y_avg += y

    y_avg /= N
    plt.plot(x, y_avg, '--', label='Average')
    plt.legend()
    plt.savefig(csvDir + 'spaghetti_plot.png')
    plt.close()

    N     = len(os.listdir(csvDir))
    y_avg = 0
    for fname in os.listdir(csvDir):
        if '.csv' not in fname:
            N -= 1
            continue
        df = pd.read_csv(csvDir + fname)
        x  = df['Time']/10
        y  = savgol_filter(df['Sliced Energy'], n, 2)
        plt.plot(x, y)

        y_avg += y

    y_avg /= N
    plt.plot(x, y_avg, '--', label='Average')
    plt.legend()
    plt.savefig(csvDir + 'spaghetti_plot_sliced.png')
    plt.close()
    return
