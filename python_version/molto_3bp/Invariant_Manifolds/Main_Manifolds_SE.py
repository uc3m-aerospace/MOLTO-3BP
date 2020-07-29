def Main_Manifolds_SE(Data):
###########################################################################
##################        SUN-EARTH MANIFOLDS       #######################
###########################################################################
# This code is devoted to compute the invariant manifolds of the Sun-Earth
# It is divided in the following sections:

### 0. Initialize General variables

### 1. SE system
#       1.1 Lagrange Points
#       1.2 Computation of periodic orbits around Libration (Lagrange) points
#           a) Obtain the initial conditions in the periodic orbit
#           b) Numerically integrate these initial conditions to obtain such periodic orbit
#           c)Apply a differential correction algorithm
#       1.3 Construct manifolds

###########################################################################
#
# Import required functions
#
    import numpy as np
    import matplotlib.pyplot as plt
    from Load_Variables import load_variables
    from IC_statemat import IC_statemat
    from Lagrange import lagrange_points, plot_lagrange_points
    from PCR3BP import PCR3BP_propagator, PCR3BP_state_derivs
    from Crossing_Det import x_crossing, y_crossing, earth_SE_crossing
    from Corrector import Corrector
    from Manifolds_tools import construct, plotm

    ## 0. Initialize General variables
    # CSpice package, Gravitational constants, Earth-Moon, and Sun-Earth constants
    print('Sun - Earth Manifolds Test Cases')

    # Define parameters (input variable to compute manifolds)
    print('Loading Variables and kernels...')
    [params_EM, params_SE] = load_variables()

    # These parameters are called by PCR3BP_state_derivs
    params = (params_SE['mu1'], params_SE['mu2'])

    ###########################################################################
    ## 1. SUN-EARTH SYSTEM
    # (Construct periodic SE orbits function)

    ## 1.1 Lagrange Points
    guess_SE = np.array([0.9, 1.01, -1])
    pos_SE   = lagrange_points(params_SE['mu2'], guess_SE) # Returns, x and y coordinates of L points

    # Plot Lagrange points
    plot_lagrange_points(params_SE['mu1'], params_SE['mu2'], pos_SE)
    # Select L2 point
    Lpoint = Data['LP'] -1
    # xe = position in x of L2
    xe = pos_SE[Lpoint, 0]

    ## 1.2 Computation of periodic orbits around Libration (Lagrange) points
    # a)Obtain the initial conditions in the periodic orbit
    # b)Numerically integrate these initial conditions to obtain such periodic orbit
    # c)Apply a differential correction algorithm

    # a)Initial condition state matrix
    mubar = params_SE['mu2']*abs(xe-1+params_SE['mu2'])**(-3) +\
        (1-params_SE['mu2'])*abs(xe+params_SE['mu2'])**(-3)

    a     = 2*mubar + 1
    b     = mubar - 1

    # Initial position in x axis: distance = x0*L_EM to L2 => Ax
    Ax = -1e-4  # = x0
                # This value says that the initial position is on the x axis,
                # at distance x0 from the libration point

    [x0, y0, vx0, vy0, eigvec, eigval, inv_phi_0] = IC_statemat(Ax, a, b)

    ## b) PCR3BP propagator: This function propagates the state of a S/C according to the PBR4BP.
    # The IC are propagated to obtain an estimation of the periodic orbit
    # Running Variables
    prnt_out_dt = Data['prnt_out_dt']   # print time period
    et0         = 0                     # t0: initial time
    deltat      = 3.5                   # time span: tf = t0 + deltat

    # stop_fun can be chosen to stop at the crossing with x or at a certain time.
    stop_fun   = x_crossing
    stop_fun.x = True

    S0 = np.array([xe+x0, y0, vx0, vy0])
    if np.imag(S0).all() == 0:
        S0 = np.real(S0)

    [SF, etf, states1_IG, times1_IG] = PCR3BP_propagator(S0, params,
        et0, deltat, prnt_out_dt, stop_fun)

    T = 2*times1_IG[-1]

    if Data['flags'][0]:

        ## c)Corrector Autonomous: Apply a differential correction algorithm using:
        # Initial conditions
        # Ax: Initial position respect L2 in x axis
        # xe: Position of L2
        # S0: Inital state vector

        # Orbital period:
        T0 = T

        # Target amplitude (target distance from xe)
        Ax_tgt      = Data['Ax_tgt']
        Ax_tgt_mid  = Ax_tgt

        # Tolerances for convergence
        Itmax   = 200
        TolRel  = 1e-10
        TolAbs  = 1e-10
        dh      = 1e-6
        Ind_Fix = 0
        Tol     = 1e-8

        dT0     = 0.1
        T0_old  = T

        exit = 0

        while abs(Ax) < abs(Ax_tgt) or exit == 0:

            print('Ax = %10.5e' % Ax)

            [X0, T0, Error, Floquet] = Corrector(PCR3BP_state_derivs, S0, params,
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
                    exit = 1

            [SF, etf, states_po, times_po] = PCR3BP_propagator (X0, params, et0, T0,
                prnt_out_dt, stop_fun)
            X0_old = X0
            T0_old = T0
            SF_old = SF
            S0 = SF
            T0 = T0 + dT0

        print('Ax = %10.5e' % Ax)
        stop_fun = 'None'
        T_po = T0_old
        prnt_out_dt = Data['prnt_out_dt']
        [SF, etf, states_po, times_po] = PCR3BP_propagator (X0_old, params, et0, T0_old,
            prnt_out_dt, stop_fun)

        # Plot corrected orbit
        plt.plot(states_po[0], states_po[1])
        plt.plot(pos_SE[1,0], pos_SE[1,1], '*')
        plt.show()

        # Save variables
        fid = open('sample.txt', 'w')
        fid.write('# Stored Data from the Corrector algorithm in Manifolds\n')
        fid.write('Orbital period: %.20f\n' % T_po)
        fid.write('# x \t\t\t\t\t y \t\t\t\t\t vx \t\t\t\t\t vy\n')
        for i in range(len(states_po[0])):
            fid.write('  %.20f %.20f %.20f %.20f\n' % (states_po[0, i],
                states_po[1, i], states_po[2, i], states_po[3, i]))
        fid.close()

    if Data['flags'][1]:
        if not Data['flags'][0]:
            fid  = open(Data['IC'], 'r')
            info = fid.readlines()
            states_po = np.array([[], [], [], []])
            for line in info:
                if line.split()[0] != '#':
                    if line.split()[:2] == ['Orbital', 'period:']:
                        T_po = np.array(float(line.split()[2]))
                    else:
                        linenum = []
                        for word in line.split():
                            linenum.append(float(word))
                        temp = np.array(linenum)
                        states_po = np.append(states_po, temp[:, np.newaxis],
                            axis = 1)
            fid.close()
            times_po = np.linspace(et0, deltat, int(deltat/prnt_out_dt) +1)



        ## 1.3 SE manifolds
        npoints = Data['npoints'] # Number of iterations = npoints*2

        stop_fun = earth_SE_crossing

        [states_s, times_s, SF_s, states_u, times_u, SF_u] = construct(
            params, T_po, states_po, times_po, eigvec, eigval,
            inv_phi_0, prnt_out_dt, npoints, stop_fun)

        plotm(params_SE['mu1'], params_SE['mu2'], pos_SE, states_po, states_s,
            SF_s, states_u)
