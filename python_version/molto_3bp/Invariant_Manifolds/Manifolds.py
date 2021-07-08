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
    from Sequential import queryFunc, queryFunctext
    from Load_Variables import load_variables
    from PCR3BP import PCR3BP_propagator, PCR3BP_state_derivs
    from Crossing_Det import poinc_crossing
    from Manifolds_tools import construct, plotm, fourierTest

    ###########################################################################

    print('\nManifolds Propagator Software\n')

    ## Initialize variables
    Data = load_variables(Data)

    ## 1. Orbit family
    if Data['input_seq']:
        if not Data['text']:
            Sequence = queryFunc()
            print('\n')
            print(Sequence)
        else:
            Sequence = queryFunctext(Data)
            print('\n')
            print(Sequence)

        from Lyapunov_Orbits.Lyapunov import Lyapunov
        from Halo_Orbits.Halo_Main import Halo_Main

        Data['flags'] = [1, 1, 1]
        Data['tf']    = 4
        Data['nmax']  = 50
        Data['tol']   = 1e-15

        orbitDef = []

        for i in range(1, Sequence['it'] +1):
            print('\nEvaluating Orbit ' + str(i))
            Data['type']   = Sequence['type' + str(i)]

            if Data['type'] == 'LY':
                Data['Ax_tgt'] = Sequence['Ax' + str(i)]
                Data['LP']     = Sequence['LP' + str(i)]
                orbitDef.append(Lyapunov(Data))
            else:
                Data['Az']  = Sequence['Az' + str(i)]
                Data['phi']  = Sequence['phi' + str(i)]
                Data['m'] = Sequence['m' + str(i)]
                Data['LP']  = Sequence['LP' + str(i)]
                orbitDef.append(Halo_Main(Data))

    else:
        if Data['type'] == 'LY':
            from Lyapunov_Orbits.Lyapunov import Lyapunov
            orbitDef = Lyapunov(Data)

        elif Data['type'] == 'HL':
            Data['flags'] = [1, 1, 1]
            Data['tf']    = 4
            Data['nmax']  = 50
            Data['tol']   = 1e-15
            from Halo_Orbits.Halo_Main import Halo_Main
            orbitDef = Halo_Main(Data)

        else:
            raise Exception('Manifolds_Main:typeError.'+\
                '    The type selected is not valid [\'LY\'][\'HL\']!')

    ## 2.1 Construct manifolds
    npoints = Data['npoints'] # Number of iterations = npoints*2

    stop_fun = poinc_crossing

    if Data['input_seq']:
        for j in range(1, Sequence['it'] +1):

            states_s = []
            SF_s     = np.array([]).reshape(-1, 4)
            states_u = []
            SF_u     = np.array([]).reshape(-1, 4)

            if 'ang' + str(j) not in Sequence:
                break

            ang = Sequence['ang' + str(j)] % 360

            if Sequence['LP' + str(j)] -1:
                angmin = min(np.arctan2(orbitDef[j-1][0][1], orbitDef[j-1][0][0]\
                    - Data['params'][0]))
                if ang*np.pi/180 > angmin +2*np.pi or ang*np.pi/180 < -angmin:
                    print('This angle is not valid (intersecting initial orbit)!')
                    if abs(ang*np.pi/180 - (angmin+2*np.pi)) < abs(ang*np.pi/180 + angmin):
                        ang = 4/3*angmin*180/np.pi +360
                        print('Reassigning to %3.2f...' % ang)
                    else:
                        ang = -4/3*angmin*180/np.pi
                        print('Reassigning to %3.2f...' % ang)
            else:
                angmin = min(np.arctan2(orbitDef[j-1][0][1], Data['params'][0]\
                    - orbitDef[j-1][0][0]))\
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

            print('Unstable manifolds...')
            [states_s_e, times_s_e, SF_s_e, states_u, times_u, SF_u] = construct(
                Data['params'], orbitDef[j-1], Data['prnt_out_dt'], npoints, 1,
                -1, stop_fun, ang, Data['pos'][Sequence['LP' + str(j)] -1][0])

            if 'Ax' + str(j +1) in Sequence or 'Az' + str(j +1) in Sequence:
                print('Stable manifolds...')
                [states_s, times_s, SF_s, states_u_e, times_u_e, SF_u_e] = construct(
                    Data['params'], orbitDef[j], Data['prnt_out_dt'], npoints, 1,
                    1, stop_fun, ang, Data['pos'][Sequence['LP' + str(j+1)] -1][0])

                print('\nPlotting manifolds\n')

                ## 2.3 Plot manifolds

                pos = np.append(Data['pos'][Sequence['LP' + str(j)] -1],
                    Data['pos'][Sequence['LP' + str(j+1)] -1])

                plotm(Data['params'][0], Data['params'][1], pos, orbitDef[j-1:j+1],
                    states_s, SF_s, states_u, SF_u, ang, angmin)

            else:

                print('\nPlotting manifolds\n')

                ## 2.3 Plot manifolds
                plotm(Data['params'][0], Data['params'][1],
                    Data['pos'][Sequence['LP' + str(j)] -1], orbitDef[j-1],
                    states_s, SF_s, states_u, SF_u, ang, angmin)

    else:
        for i in Data['poincSec']:

            ang = i % 360
            if Data['LP'] -1:
                angmin = min(np.arctan2(orbitDef[0][1], orbitDef[0][0] - Data['params'][0]))
                if ang*np.pi/180 > angmin +2*np.pi or ang*np.pi/180 < -angmin:
                    print('This angle is not valid (intersecting initial orbit)!')
                    if abs(ang*np.pi/180 - (angmin+2*np.pi)) < abs(ang*np.pi/180 + angmin):
                        ang = 4/3*angmin*180/np.pi +360
                        print('Reassigning to %3.2f...' % ang)
                    else:
                        ang = -4/3*angmin*180/np.pi
                        print('Reassigning to %3.2f...' % ang)
            else:
                angmin = min(np.arctan2(orbitDef[0][1], Data['params'][0] - orbitDef[0][0]))\
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
                Data['params'], orbitDef, Data['prnt_out_dt'], npoints, Data['d'],
                Data['branch'], stop_fun, ang, Data['pos'][Data['LP'] -1][0])

            print('\nPost-processing Data...\n')

            ## 2.2 Fourier analysis
            fourierTest(Data['params'][0], Data['params'][1], Data['pos'][Data['LP'] -1],
                states_s, states_u, ang, Data)

            print('\nPlotting manifolds\n')

            ## 2.3 Plot manifolds
            plotm(Data['params'][0], Data['params'][1], Data['pos'][Data['LP'] -1],
                orbitDef, states_s, SF_s, states_u, SF_u, ang, angmin)
