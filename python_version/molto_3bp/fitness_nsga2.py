def fitness_nsga2(x, params):
    # Import functions required
    import math
    from PCRFBP_propagator import PCRFBP_propagator
    ## Identify decision variables
    alpha = x[0]
    beta  = x[1]
    ti    = x[2]
    delta = x[3]
    #
    r0 = params['r0']
    rf = params['rf']
    mu = params['mu']
    DU = params['DU']
    VU = params['VU']
    #
    prnt_out_dt = 0.001
    stop_fun = 'none'
    #
    ## Initial state
    x0    = r0*math.cos(alpha)-mu
    y0    = r0*math.sin(alpha)
    v0    = beta*math.sqrt((1.0-mu)/r0)
    x0dot = -(v0-r0)*math.sin(alpha)
    y0dot = (v0-r0)*math.cos(alpha)

    S0 = [x0, y0, x0dot, y0dot]

    ## Propagate PCRFBP
    [SF, etf, states, times] = PCRFBP_propagator (S0, params,
        ti, delta, prnt_out_dt, stop_fun)

    xf = SF[0]
    yf = SF[1]
    xfdot = SF[2]
    yfdot = SF[3]
    alpha = math.atan2(yf,xf+mu-1)
    Vm = math.sqrt(mu/rf)

    ## Evaluate objective functions
    f = [0, 0]
    f[0] = math.fabs(math.sqrt((xf+mu-1.0)**2+yf**2)-rf)*DU/1000.0 # delta r (km - dimensional units)
    f[1] = math.fabs(math.sqrt((-Vm*math.sin(alpha)-(xfdot-yf))**2 \
    + (Vm*math.cos(alpha)-(yfdot+xf+mu-1.0))**2))*VU/1000.0; # velocity vector error (km/s - dimensional units)

    cons = 0
    return f, cons
