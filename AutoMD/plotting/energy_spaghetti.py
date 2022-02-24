from ..config import *
from ast import literal_eval

def energy_spaghetti(csvDir, n=51):
    def prep_data(df, property):
        tmp = {}
        for i in range(len(df)):
            e = literal_eval(df[property][i])
            for j in range(len(e)):
                name = 'mol' + str(j+1)
                if name not in tmp:
                    tmp[name] = []
                tmp[name].append(e[j])

        df2 = pd.DataFrame(tmp)
        df2 = pd.concat([df2, df['Time']], axis=1)
        return df2

    def loop(property, ext):
        N     = len(os.listdir(csvDir))
        y_avg = 0
        for fname in os.listdir(csvDir):
            if '.csv' not in fname:
                N -= 1
                continue
            df = pd.read_csv(csvDir + fname)
            df = prep_data(df, property)
            x  = df['Time']/10
            for key in df:
                if 'Time' in key:
                    continue
                y  = savgol_filter(df[key], n, 2)
                plt.plot(x, y)

            plt.savefig(csvDir + fname.split('_')[0] + ext)
            plt.close()
        return
    
    loop('All Vibra Energy', 'vibra_spaghetti.png')
    loop('All Rotat Energy', 'rotat_spaghetti.png')
    loop('All Trans Energy', 'trans_spaghetti.png')
    return
