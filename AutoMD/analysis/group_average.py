from ..config import *

def group_average(h5):
    hdf = pd.HDFStore(h5)
    tmp = {}

    # i is [groupName, groupsInGroup, membersInGroup]
    for i in hdf.walk():
        for j in i[-1]:

            # key becomes: groupName - last char + 'Avg'
            key = f'{i[0][:-1]}Avg'

            # Here we just fill dictionary with
            # strings for getting memebers
            try:
                tmp[key].append(f'{i[0]}/{j}') # ie. groupName/member
            except:
                tmp[key] = [f'{i[0]}/{j}']

    # Now we run through the dictionary
    # to get members and average the
    # DataFrame into a zero-ed out one
    for key in tmp.keys():
        avg = hdf.get(tmp[key][0]) * 0

        for i in tmp[key]:
            df   = hdf.get(i)
            avg += df

        avg /= len(tmp[key])

        hdf.put(key, df)

    hdf.close()
    return
