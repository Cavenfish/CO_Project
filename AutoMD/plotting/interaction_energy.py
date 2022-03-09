from ..config import *

def interaction_energy(csvDir):
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

    x = avg['Time']
    y = avg['Total Energy'] - avg['Sliced Energy']
    plt.plot(x,y)
    plt.show()

    return