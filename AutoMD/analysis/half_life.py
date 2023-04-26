from ..config import *

#Single exponential decay
def f(x, E, tau):
    return (E * np.exp(-x / tau))

#Two channel exponential decay
def g(x, E, tau, E2, tau2):
    return f(x, E, tau) + f(x, E2, tau2)

#The two channel expoential decay fits the
#decay curves best. The assumption here is:
#-> One channel is for vib-vib transfer
#-> Other is for vib-lattice transfer

def get_fit(df, p0):
    dp0 = np.array([0.01, 500, -0.01, 30])
    x   = df['Time']
    y   = df['Sliced Energy']

    try:
        bounds     = ([0,0,0,0],[1,1e10,1,1e4])
        popt, pcov = curve_fit(g, x, y, p0, bounds=bounds)
    except RuntimeError:
        p0 += dp0
        popt, pcov = get_fit(df, p0)

    return popt, pcov

def half_life(h5):
    hdf = pd.HDFStore(h5)

    for i in hdf.walk():

        # Prevent entering in unwanted groups
        if (len(i[-1]) < 20) or ('VACF' in i[0]):
            continue

        tmp = {'uid':[], 'E1':[], 'tau1':[], 'E2':[], 'tau2':[]}
        for j in i[-1]:
            if not j.isdigit():
                continue

            df = hdf.get(f'{i[0]}/{j}')
            p0 = np.array([0.01, 1000, 0.4, 1])

            popt, pcov = get_fit(df, p0)

            tmp['uid' ].append(j)
            tmp['E1'  ].append(popt[0])
            tmp['tau1'].append(popt[1] * np.log(2))
            tmp['E2'  ].append(popt[2])
            tmp['tau2'].append(popt[3] * np.log(2))

        df = pd.DataFrame(tmp)
        hdf.put(f'{i[0]}/tau', df)
    hdf.close()
    return
