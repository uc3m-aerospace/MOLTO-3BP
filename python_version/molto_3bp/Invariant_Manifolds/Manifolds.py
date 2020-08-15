def Manifolds(Data):
###########################################################################
#######################        MANIFOLDS       ############################
###########################################################################
# This code is devoted to compute the invariant manifolds of any input orbit
# It is divided in the following sections:

### 0. Initialize common variables

### 1. Initialize Orbit family

### 2. Orbit type (2D or 3D)
#       2.1 Construct manifolds
#       2.2 Fourier analysis
#       2.3 Plot manifolds

###########################################################################
#
# Import required functions
#
    import numpy as np
    from Load_Variables import load_variables
    from PCR3BP import PCR3BP_propagator, PCR3BP_state_derivs
    from Crossing_Det import poinc_crossing
    from Manifolds_tools import construct, plotm, fourierTest

    ###########################################################################

    print('\nManifolds Propagator Software\n')

    ## Initialize variables
    Data = load_variables(Data)

    ## 1. Orbit family
    if Data['type'] == 'LY':
        from Lyapunov_Orbits.Lyapunov import Lyapunov
        (states_po, times_po, T_po, eigvec) = Lyapunov(Data)

    elif Data['type'] == 'HL':
        Data['flags'] = [1, 1, 1]
        Data['tf']    = 4
        Data['nmax']  = 50
        Data['tol']   = 1e-15
        from Halo_Orbits.Halo_Main import Halo_Main
        (states_po, times_po, T_po, eigvec) = Halo_Main(Data)

    else:
        raise Exception('Manifolds_Main:typeError.'+\
            '    The type selected is not valid [\'LY\'][\'HL\']!')

    ## 2.1 Construct manifolds
    npoints = Data['npoints'] # Number of iterations = npoints*2

    stop_fun = poinc_crossing

    for i in Data['poincSec']:

        ang = i % 360
        if Data['LP'] -1:
            angmin = min(np.arctan2(states_po[1], states_po[0] - Data['params'][0]))
            if ang*np.pi/180 > angmin +2*np.pi or ang*np.pi/180 < -angmin:
                print('This angle is not valid (intersecting initial orbit)!')
                if abs(ang*np.pi/180 - (angmin+2*np.pi)) < abs(ang*np.pi/180 + angmin):
                    ang = 4/3*angmin*180/np.pi +360
                    print('Reassigning to %3.2f...' % ang)
                else:
                    ang = -4/3*angmin*180/np.pi
                    print('Reassigning to %3.2f...' % ang)
        else:
            angmin = min(np.arctan2(states_po[1], Data['params'][0] - states_po[0]))\
                + np.pi
            if ang*np.pi/180 > angmin and ang*np.pi/180 < -angmin +2*np.pi:
                print('This angle is not valid (intersecting initial orbit)!')
                if abs(ang*np.pi/180 - (-angmin+2*np.pi)) < abs(ang*np.pi/180 - angmin):
                    ang = -4/3*angmin*180/np.pi + 420
                    print('Reassigning to %3.2f...' % ang)
                else:
                    ang = 4/3*angmin*180/np.pi - 60
                    print('Reassigning to %3.2f...' % ang)

        print('\nConstructing Manifolds...\n')
        print('Poincaré section angle = %3.1f' % ang)

        [states_s, times_s, SF_s, states_u, times_u, SF_u] = construct(
            Data['params'], T_po, states_po, times_po, eigvec,
            Data['prnt_out_dt'], npoints, Data['d'], Data['branch'],
            stop_fun, ang, Data['pos'][Data['LP'] -1][0])

        print('\nPost-processing Data...\n')

        ## 2.2 Fourier analysis
        fourierTest(Data['params'][0], Data['params'][1], Data['pos'][Data['LP'] -1],
            states_s, states_u, ang, Data)

        print('\nPlotting manifolds\n')

        ## 2.3 Plot manifolds
        plotm(Data['params'][0], Data['params'][1], Data['pos'][Data['LP'] -1],
            states_po, states_s, SF_s, states_u, SF_u, ang, angmin)
