def x_crossing(et, state, mu1, mu2):
    return state[1]

def y_crossing(et, state, mu1, mu2):
    return state[0]

def poinc_crossing(et, state, mu1, mu2, ang, L):

    import numpy as np

    [x, y] = state[:2]

    d = abs((L - mu1)*np.sin(ang*np.pi/180)*x\
        - (L - mu1)*np.cos(ang*np.pi/180)*y\
        - (L - mu1)*np.sin(ang*np.pi/180)*mu1)\
        /(np.sqrt(((L - mu1)*np.sin(ang*np.pi/180))**2\
        + ((L - mu1)*np.cos(ang*np.pi/180))**2))

    if np.log10(mu2) < -5:
        tol = 1e-5
    else:
        tol = 5e-4

    if L > mu1:
        cond = x < L
    else:
        cond = x > L

    if d < tol and np.sign(np.cos(ang*np.pi/180)) == np.sign(x-mu1) and\
        np.sign(np.sin(ang*np.pi/180)) == np.sign(np.round(y, decimals = 4)) and\
        cond:
        return 0
    elif abs(y) > abs(L-mu1):
        return 0
    else:
        return 1
