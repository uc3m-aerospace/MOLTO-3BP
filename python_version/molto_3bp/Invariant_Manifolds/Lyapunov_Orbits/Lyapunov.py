def Lyapunov(Data):
###########################################################################
####################        Lyapunov Orbits       #########################
###########################################################################
# This code is devoted to compute members of the lyapunov family of orbits
# It is divided in the following sections:

### 0. Initialize General variables

### 1. SE or EM system
#       1.1 Computation of periodic orbits around Libration (Lagrange) points
#           a) Obtain the initial conditions in the periodic orbit
#           b) Numerically integrate these initial conditions to obtain such
#               periodic orbit
#           c)Apply a differential correction algorithm

###########################################################################
#
# Import required functions
#
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from .IC_statemat import IC_statemat
    from .PCR3BP import PCR3BP_propagator, PCR3BP_state_derivs
    from .Corrector import Corrector

    print('Lyapunov Generator\n')

    ###########################################################################
    ## 1. SE or EM SYSTEM
    # (Construct periodic orbits function)

    # xe = position in x of L
    # Select L point
    Lpoint = Data['LP'] -1

    xe = Data['pos'][Lpoint, 0]

    ## 1.1 Computation of periodic orbits around Libration (Lagrange) points
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
        Ax = -1e-3*np.sign(xe - Data['params'][0])
                # This value says that the initial position is on the x axis,
                # at distance x0 from the libration point

    [x0, y0, vx0, vy0] = IC_statemat(Ax, a, b)

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

        [X0, T0, Error, eigvec] = Corrector(PCR3BP_state_derivs, S0, Data['params'],
            T0, Itmax, Tol, TolRel, TolAbs, dh, Ind_Fix)

        # X0 are the corrected IC and T0 is the corrected period
        # T0            # Corrected period
        Ax = X0[0]-xe   # Distance to L point corrected

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
    ax.plot(Data['pos'][Lpoint,0], Data['pos'][Lpoint,1], 'k*')
    plt.title('Corrected Orbit in the Synodic ref. frame')
    plt.xlabel('Distance from baricenter in nondimensional coordinates')
    plt.ylabel('Distance perpendicular to the line of primaries')
    ax.ticklabel_format(style = 'sci', scilimits = (-3, 3), useMathText = True)
    plt.show()

    # Save variables
    fid = open("Lyapunov_Orbits" + os.sep + 'sample.txt', 'w')
    fid.write('# Stored Data from the Corrector algorithm in Lyapunov\n')
    fid.write('Orbital period: %.20f\n' % T_po)
    fid.write('# x0 \t\t\t\t\t y0 \t\t\t\t\t vx0 \t\t\t\t\t vy0\n')
    fid.write('  %.20f %.20f %.20f %.20f\n' % (states_po[0, 0],
        states_po[1, 0], states_po[2, 0], states_po[3, 0]))
    fid.close()

    return (states_po, T_po, eigvec)
