def poinc_crossing(et, state, mu1, mu2, ang, L):

    import numpy as np

    [x, y] = state[:2]

    d = ((L - mu1)*np.sin(ang*np.pi/180)*(x-mu1)\
        - (L - mu1)*np.cos(ang*np.pi/180)*y)

    xprime = np.cos(ang*np.pi/180)*(x-mu1) + np.sin(ang*np.pi/180)*y

    cond = xprime > 0

    if abs(y) > abs(L-mu1):
        return 0
    elif cond:
        return d
    else:
        if ang - 180 > 0:
            return -mu1
        else:
            return mu1
