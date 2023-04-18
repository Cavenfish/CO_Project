from ..config import *

def group_average(h5):
    hdf = pd.HDFStore(h5) 

    tmp = {}
    for i in hdf.walk():
        for j in i[-1]:
            key = f'{i[0][:-1]}Avg'
            
            try:
                tmp[key].append(f'{i[0]}/{j}')
            except:
                tmp[key] = [f'{i[0]}/{j}']


    for key in tmp.keys():
        avg = hdf.get(tmp[key][0]) * 0
        
        for i in tmp[key]:
            df   = hdf.get(i)
            avg += df
        
        avg /= len(tmp[key])
        
        hdf.put(key, df)

    hdf.close()
    return
