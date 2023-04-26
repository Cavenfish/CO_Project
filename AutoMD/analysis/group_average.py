from ..config import *

def leaves(hdf, s):
    for i in range(2):
        node = hdf.get_node(s + f'{i}')

        try:
            members = node.__members__
        except:
            continue

        for j in members:
            if 'tau' in j:
                continue
            yield f'{s}{i}/{j}'

def group_average(h5):
    hdf = pd.HDFStore(h5)

    # i is [groupName, groupsInGroup, membersInGroup]
    for i in hdf.walk():
        if ('co0' not in i[1]) and ('13c18o0' not in i[1]):
            continue

        for clu in ['co', '13c18o']:
            n   = 0
            key = f'{i[0]}/{clu}'
            print(key)
            for leaf in leaves(hdf, key):
                df = hdf.get(leaf)

                if n == 0:
                    avg = df * 0

                avg += df
                n   += 1

            avg /= n
            hdf.put(key+'Avg', avg)

    hdf.close()
    return
