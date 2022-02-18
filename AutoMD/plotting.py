from .config import *

def spaghetti_plot(csvDir, n=5001):
    
    n     = len(os.listdir(csvDir))
    y_avg = 0 
    for fname in os.listdir(csvDir):
        df = pd.read_csv(fname)
        x  = df['Time']/1000
        y  = savgol_filter(df['Total Energy'], n, 2)
        plt.plot(x, y)
        
        y_avg += y

    y_avg /= n
    plt.plot(x, y_avg, '--', label='Average')
    return

