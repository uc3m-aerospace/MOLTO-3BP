def x_crossing(et, state, mu1, mu2):
    return state[1]

def y_crossing(et, state, mu1, mu2):
    if et > 0.:
        return state[0]
    else:
        return 1
