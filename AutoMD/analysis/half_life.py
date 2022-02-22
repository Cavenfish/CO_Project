from ..config import *

def half_life(directory):
    tauList = []
    for filename in os.listdir(directory):
        df     = pd.read_csv(filename)
        E, tau = half_life(df, filename.replace('.csv', '.png'))
        tauList.append(tau)

    tau_avg = np.average(tauList)

    return


