def Manifolds(Data):
###########################################################################
#######################        MANIFOLDS       ############################
###########################################################################
# This code is devoted to compute the invariant manifolds of any input orbit
# It is divided in the following sections:

### 0. Initialize Orbit family

### 1. Orbit type (2D or 3D)
#       1.1 Construct manifolds
#       1.2 Fourier analysis
#       1.3 Plot manifolds

###########################################################################
#
# Import required functions
#
    import numpy as np
    from PCR3BP import PCR3BP_propagator, PCR3BP_state_derivs
    from Crossing_Det import poinc_crossing
    from Manifolds_tools import construct, plotm, fourierTest

    ###########################################################################

    ## 0. Orbit family
    if Data['type'] == 'LY':
        from Lyapunov_Orbits.Lyapunov import Lyapunov
        (Data, states_po, times_po, T_po, eigvec, eigval, inv_phi_0, pos) =\
            Lyapunov(Data)
        [xL, yL] = pos[Data['LP'] -1]

    elif Data['type'] == 'HL':
        Data['flags'] = [1, 1, 1]
        Data['tf']    = 4
        Data['nmax']  = 50
        Data['tol']   = 1e-15
        from Halo_Orbits.Halo_Main import Halo_Main
        (Data, states_po, times_po, T_po, eigvec, eigval, inv_phi_0, xL) =\
            Halo_Main(Data)
        xL = xL[0]
        yL = 0
    else:
        raise Exception('Manifolds_Main:typeError.'+\
            '    The type selected is not valid [\'LY\'][\'HL\']!')

    print('\nConstructing Manifolds...\n')

    ## 1.1 Construct manifolds
    npoints = Data['npoints'] # Number of iterations = npoints*2

    stop_fun = poinc_crossing


    ang = Data['poincSec'] % 360
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

    [states_s, times_s, SF_s, states_u, times_u, SF_u] = construct(
        Data['params'], T_po, states_po, times_po, eigvec, eigval,
        inv_phi_0, Data['prnt_out_dt'], npoints, Data['d'], stop_fun, ang,
        xL)

    print('\nPost-processing Data...\n')

    ## 1.2 Fourier analysis
    fourierTest(Data['params'][0], Data['params'][1], [xL, yL],
        states_s, states_u, ang, Data)

    print('\nPlotting manifolds\n')

    ## 1.3 Plot manifolds
    plotm(Data['params'][0], Data['params'][1], [xL, yL], states_po,
        states_s, SF_s, states_u, SF_u, ang, angmin)
