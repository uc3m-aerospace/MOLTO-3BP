def Lyapunov(Data):
###########################################################################
####################        Lyapunov Orbits       #########################
###########################################################################
# This code is devoted to compute members of the lyapunov family of orbits
# It is divided in the following sections:

### 0. Initialize General variables

### 1. SE or EM system
#       1.1 Lagrange Points
#       1.2 Computation of periodic orbits around Libration (Lagrange) points
#           a) Obtain the initial conditions in the periodic orbit
#           b) Numerically integrate these initial conditions to obtain such
#               periodic orbit
#           c)Apply a differential correction algorithm

###########################################################################
#
# Import required functions
#
    import numpy as np
    import matplotlib.pyplot as plt
    from .Load_Variables import load_variables
    from .IC_statemat import IC_statemat
    from .Lagrange import lagrange_points, plot_lagrange_points
    from .PCR3BP import PCR3BP_propagator, PCR3BP_state_derivs
    from .Corrector import Corrector

    ## 0. Initialize General variables
    # CSpice package, Gravitational constants, Earth-Moon, and Sun-Earth constants
    print('SE & EM Manifolds Test Cases Software')
    if Data['mode'] not in ['SE', 'EM']:
        raise Exception('Manifolds_Main:modeError.'+\
            '    The mode selected is not valid [\'SE\'][\'ME\']!')

    # Define parameters (input variable to compute manifolds)
    print('Loading Variables and kernels...')
    [params_EM, params_SE] = load_variables()

    # These parameters are called by PCR3BP_state_derivs
    if Data['mode'] == 'SE':
        Data['params'] = (params_SE['mu1'], params_SE['mu2'])
    else:
        Data['params'] = (params_EM['mu1'], params_EM['mu2'])

    ###########################################################################
    ## 1. SE or EM SYSTEM
    # (Construct periodic orbits function)

    ## 1.1 Lagrange Points
    guess = np.array([0.9, 1.01, -1])
    pos   = lagrange_points(Data['params'][1], guess) # Returns, x and y coordinates of L points

    # Plot Lagrange points
    plot_lagrange_points(Data['params'][0], Data['params'][1], pos)
    # Select L2 point
    Lpoint = Data['LP'] -1
    # xe = position in x of L2
    xe = pos[Lpoint, 0]

    ## 1.2 Computation of periodic orbits around Libration (Lagrange) points
    # a)Obtain the initial conditions in the periodic orbit
    # b)Numerically integrate these initial conditions to obtain such periodic orbit
    # c)Apply a differential correction algorithm

    # a)Initial condition state matrix
    mubar = Data['params'][1]*abs(xe-1+Data['params'][1])**(-3) +\
        (1-Data['params'][1])*abs(xe+Data['params'][1])**(-3)

    a     = 2*mubar + 1
    b     = mubar - 1

    # Initial position in x axis: distance = x0*L_EM to L2 => Ax
    if Data['mode'] == 'SE':
        Ax = -1e-4*np.sign(xe - Data['params'][0])  # Non dimensional distance
    else:
        if Lpoint:
            Ax = -2e-3
        else:
            Ax = 1e-3
                # This value says that the initial position is on the x axis,
                # at distance x0 from the libration point

    [x0, y0, vx0, vy0, eigvec, eigval, inv_phi_0] = IC_statemat(Ax, a, b)

    ## b) PCR3BP propagator: This function propagates the state of a S/C
    # according to the PBR3BP.
    # The IC are propagated to obtain an estimation of the periodic orbit
    # Running Variables
    prnt_out_dt = Data['prnt_out_dt']   # print time period
    et0         = 0                     # t0: initial time
    if Data['mode'] == 'SE':
        deltat = 3.5    # time span: tf = t0 + deltat
    else:
        deltat = 10

    # stop_fun can be chosen to stop at the crossing with x axis
    stop_fun   = lambda et, state, mu1, mu2: state[1]
    stop_fun.x = True

    S0 = np.array([xe+x0, y0, vx0, vy0])
    if np.imag(S0).all() == 0:
        S0 = np.real(S0)

    [SF, etf, states1_IG, times1_IG] = PCR3BP_propagator(S0, et0, deltat,
        prnt_out_dt, stop_fun, Data['params'][0], Data['params'][1])

    T = 2*times1_IG[-1]

    ## c)Corrector Autonomous: Apply a differential correction algorithm using:
    # Initial conditions
    # Ax: Initial position respect L2 in x axis
    # xe: Position of L2
    # S0: Inital state vector

    # Orbital period:
    T0 = T

    # Target amplitude (target distance from xe)
    Ax_tgt      = Data['Ax_tgt']
    Ax_tgt_mid  = np.copy(Ax_tgt)

    # Tolerances for convergence
    Itmax   = 200
    TolRel  = 1e-10
    TolAbs  = 5e-11
    dh      = 1e-6
    Ind_Fix = 0
    Tol     = 1e-8

    dT0     = 0.1
    T0_old  = T

    stop_fun.x = False
    stop_fun.direction = 0

    while abs(Ax) < abs(Ax_tgt):

        print('Ax = %10.5e & Ax target = %10.5e' % (abs(Ax), abs(Ax_tgt)))

        [X0, T0, Error, Floquet] = Corrector(PCR3BP_state_derivs, S0, Data['params'],
            T0, Itmax, Tol, TolRel, TolAbs, dh, Ind_Fix)

        # X0 are the corrected IC and T0 is the corrected period
        # T0            # Corrected period
        Ax = X0[0]-xe   # Distance to L2 corrected

        if (T0 < 1 or abs(Ax) > abs(Ax_tgt_mid) or T0 > T0_old*2):
            X0  = X0_old
            T0  = T0_old
            dT0 = dT0/2
            Ax_tgt_mid = 2*Ax_tgt

            if dT0 < 1e-3 or T0 > T0_old*2 or abs(Ax) > abs(Ax_tgt):
                break

        [SF, etf, states_po, times_po] = PCR3BP_propagator (X0, et0, T0,
            prnt_out_dt*10, stop_fun, Data['params'][0], Data['params'][1])

        X0_old = X0
        T0_old = T0
        SF_old = SF
        S0 = SF
        T0 = T0 + dT0

    print('Ax = %10.5e' % Ax)
    stop_fun = 'None'
    T_po = T0_old
    prnt_out_dt = Data['prnt_out_dt']
    [SF, etf, states_po, times_po] = PCR3BP_propagator (X0_old, et0, T0_old,
        prnt_out_dt, stop_fun, Data['params'][0], Data['params'][1])

    # Plot corrected orbit
    fig, ax = plt.subplots()
    ax.plot(states_po[0], states_po[1])
    ax.plot(pos[Lpoint,0], pos[Lpoint,1], 'k*')
    plt.title('Corrected Orbit in the Synodic ref. frame')
    plt.xlabel('Distance from baricenter in nondimensional coordinates')
    plt.ylabel('Distance perpendicular to the line of primaries')
    ax.ticklabel_format(style = 'sci', scilimits = (-3, 3), useMathText = True)
    plt.show()

    # Save variables
    fid = open('sample.txt', 'w')
    fid.write('# Stored Data from the Corrector algorithm in Lyapunov\n')
    fid.write('Orbital period: %.20f\n' % T_po)
    fid.write('# x0 \t\t\t\t\t y0 \t\t\t\t\t vx0 \t\t\t\t\t vy0\n')
    fid.write('  %.20f %.20f %.20f %.20f\n' % (states_po[0, 0],
        states_po[1, 0], states_po[2, 0], states_po[3, 0]))
    fid.close()

    return (Data, states_po, times_po, T_po, eigvec, eigval, inv_phi_0, pos)
