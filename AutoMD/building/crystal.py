from ..config import *

#This code is so fucking hacky

#TODO: Make it less hacky lol

def p213(x, y, z, abc):
    pos = []
    l   = [[    x,     y,     z],
           [0.5+x, 0.5-y,    -z],
           [   -x, 0.5+y, 0.5-z],
           [0.5-x,    -y, 0.5+z],
           [    z,     x,     y],
           [0.5-z,    -x, 0.5+y],
           [0.5+z, 0.5-x,    -y],
           [   -z, 0.5+x, 0.5-y],
           [    y,     z,     x],
           [   -y, 0.5+z, 0.5-x],
           [0.5-y,    -z, 0.5+x],
           [0.5+y, 0.5-z,    -x]]

    for i in l:
        a = []
        for j in i:
            j *= abc
            if j > abc:
                j -= abc
            if j < 0:
                j += abc
            a.append(j)

        if a in pos:
            continue
        flag = False
        for k in pos:
            diff = np.array(k) - np.array(a)
            d    = np.sqrt(np.dot(diff,diff))
            if d < 1:
                flag = True
                break
        if flag:
            continue
        pos.append(a)
    
    return pos
