def Manifolds(Data):
###########################################################################
#######################        MANIFOLDS       ############################
###########################################################################
# This code is devoted to compute the invariant manifolds of the Sun-Earth
# It is divided in the following sections:

### 0. Initialize General variables

### 1. SE or EM system
#       1.1 Lagrange Points
#       1.2 Computation of periodic orbits around Libration (Lagrange) points
#           a) Obtain the initial conditions in the periodic orbit
#           b) Numerically integrate these initial conditions to obtain such
#               periodic orbit
#           c)Apply a differential correction algorithm
#       1.3 Construct manifolds
#       1.4 Plot manifolds

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
    from Crossing_Det import x_crossing, y_crossing, poinc_crossing
    from Corrector import Corrector
    from Manifolds_tools import construct, plotm

    ## 0. Initialize General variables
    # CSpice package, Gravitational constants, Earth-Moon, and Sun-Earth constants
    print('SE & EM Manifolds Test Cases Software')
    if Data['mode'] not in ['SE', 'EM', 'ALL']:
        raise Exception('Manifolds_Main:modeError.'+\
            '    The mode selected is not valid [\'SE\'][\'ME\'][\'ALL\']!')

    # Define parameters (input variable to compute manifolds)
    print('Loading Variables and kernels...')
    [params_EM, params_SE] = load_variables()

    # These parameters are called by PCR3BP_state_derivs
    if Data['mode'] == 'SE':
        params = (params_SE['mu1'], params_SE['mu2'])
    else:
        params = (params_EM['mu1'], params_EM['mu2'])

    ###########################################################################
    ## 1. SE or EM SYSTEM
    # (Construct periodic orbits function)

    ## 1.1 Lagrange Points
    guess = np.array([0.9, 1.01, -1])
    pos   = lagrange_points(params[1], guess) # Returns, x and y coordinates of L points

    # Plot Lagrange points
    plot_lagrange_points(params[0], params[1], pos)
    # Select L2 point
    Lpoint = Data['LP'] -1
    # xe = position in x of L2
    xe = pos[Lpoint, 0]

    ## 1.2 Computation of periodic orbits around Libration (Lagrange) points
    # a)Obtain the initial conditions in the periodic orbit
    # b)Numerically integrate these initial conditions to obtain such periodic orbit
    # c)Apply a differential correction algorithm

    # a)Initial condition state matrix
    mubar = params[1]*abs(xe-1+params[1])**(-3) +\
        (1-params[1])*abs(xe+params[1])**(-3)

    a     = 2*mubar + 1
    b     = mubar - 1

    # Initial position in x axis: distance = x0*L_EM to L2 => Ax
    if Data['mode'] == 'SE':
        Ax = -1e-4  # Non dimensional distance
    else:
        Ax = -2e-3
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

    # stop_fun can be chosen to stop at the crossing with x or at a certain time.
    stop_fun   = x_crossing
    stop_fun.x = True

    S0 = np.array([xe+x0, y0, vx0, vy0])
    if np.imag(S0).all() == 0:
        S0 = np.real(S0)

    [SF, etf, states1_IG, times1_IG] = PCR3BP_propagator(S0, et0, deltat,
        prnt_out_dt, stop_fun, params[0], params[1])

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
                    break

            [SF, etf, states_po, times_po] = PCR3BP_propagator (X0, et0, T0,
                prnt_out_dt*10, stop_fun, params[0], params[1])

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
            prnt_out_dt, stop_fun, params[0], params[1])

        # Plot corrected orbit
        plt.plot(states_po[0], states_po[1])
        plt.plot(pos[1,0], pos[1,1], '*')
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
            times_po = np.linspace(et0, T_po, int(T_po/prnt_out_dt) +1)

        ## 1.3 Construct manifolds
        npoints = Data['npoints'] # Number of iterations = npoints*2

        stop_fun = poinc_crossing

        [states_s, times_s, SF_s, states_u, times_u, SF_u] = construct(
            params, T_po, states_po, times_po, eigvec, eigval,
            inv_phi_0, prnt_out_dt, npoints, stop_fun, Data['poincSec'] % 360,
            pos[Lpoint, 0])

        ## 1.4 Plot manifolds
        plotm(params[0], params[1], pos, states_po, states_s,
            SF_s, states_u, SF_u, Data['poincSec'] % 360)
