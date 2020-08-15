def Halo_Num_Comp(Data):
    ##### HALO ORBITS NUMERICAL COMPUTATION #####
    #
    # Importing required functions
    #
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy import linalg
    from .intFun import DiffCorrection

    ## IC guess

    # Use Results from Halo_Generator.py or sample.txt

    #--- Initial Guess ---
    if Data['flags'][0] or Data['method'] == 'insitu':
        [x0, z0, vy0] = Data['IC']
    else:
        fid  = open(Data['IC'],'r')
        info = fid.readlines()
        IC = []
        for i in info:
            if i.split()[0] == '#':
                IC.append(i.split()[-1])
        if len(IC) == 3:
            [x0, z0, vy0] = IC
        elif len(IC) == 6:
            [x0, z0, vy0] = [IC[0], IC[2], IC[4]]
        else:
            raise Exception('Halo_Num_Comp:ICError.' +\
                '   The text file selected does not have the right format!')

    # Recall y0, vx0, vz0 = 0

    ## Differential Correction
    for i in range(1, Data['nmax']+1):

    # IC vector preparation

    # x' = f(x)    IC
        q0 = np.zeros(42) # Initial Conditions
        q0[:6] = [x0, 0, z0, 0, vy0, 0]

    # Phi' = Df(x)*Phi(t,t0)     IC
        Phi0 = np.identity(6)
        q0[6:] = Phi0.ravel()

    # Ode45 - State Transition Matrix Computation

    # Event to stop ode45, x-z plane cross (at T/2)
        def EventHalo(t, q, mu):
            if t > 0.:
                return q[1]
            else:
                return 1
        EventHalo.terminal = True
        sol = solve_ivp(DiffCorrection, [0, Data['tf']],
            q0, events = EventHalo, args = (Data['mu'],),
            atol = 1e-8,
            rtol = 1e-5)
        q = sol.y
        t = sol.t

    # Extracting solution
        xfvec  = q[:6,-1]
        xbfvec = q[:6,-2]

        Phifvec = q[6:,-1]
        Phifmat = Phifvec.reshape(6, -1) # to matrix form

    # Desired values at tf: vxf,vzf = 0   (yf=0 already with stop event)
    # Desired values             % Final values
        vxfdes = 0; vzfdes = 0; vxf = xfvec[3]; vzf = xfvec[5];

        ydot = xfvec[4]
        xdotdot = (vxf-xbfvec[3])/(t[-1]-t[-2])
        zdotdot = (vzf-xbfvec[5])/(t[-1]-t[-2])

    # Delta x
        dvx = vxfdes-vxf; dvz = vzfdes-vzf
        B = np.array([[dvx],[dvz]])
        D = np.array([[xdotdot],[zdotdot]])
        E = np.array([Phifmat[1,0], Phifmat[1,4]])

    # Check of IC

        err1 = abs(dvx); err2 = abs(dvz)

        if (err1 <= Data['tol']) and (err2 <= Data['tol']):
            break
        else:

    # Update IC --- Ax=B
            A = np.array(Phifmat[np.array([3, 3, 5, 5]),
                np.array([0, 4, 0, 4])].reshape((2,2)))
            C = A-1/ydot*D*E
            dxvec0 = linalg.solve(C,B) # Solve inverting C
            dx0    = dxvec0[0]
            dvy0   = dxvec0[1]

            x0  =  x0 + dx0
            vy0 = vy0 + dvy0

    ## Solution

    print('--- Halo Generator: Numerical Computation ---\n')
    if (err1 <= Data['tol']) and (err2 <= Data['tol']):
        print('Nº of iterations to converge: ' + str(i))
        print('\n--- Solution ---')
        print('x0  = %.20f;' % x0)
        print('y0  = 0.00;')
        print('z0  = %.20f;' % z0)
        print('vx0 = 0.00;')
        print('vy0 = %.20f;' % vy0)
        print('vz0 = 0.00;\n')
        print('--- Orbit Period ---')
        print('T = %.20f;' % (t[-1]*2))
        print('T/2 = %.20f;\n' % t[-1])
    else:
        print('The program has not converged!')
        print('err1  = ' + str(err1))
        print('err2  = ' + str(err2))
        print('\nNº of iterations done: ' + str(i))
        print('Tolerance: ' + str(Data['tol']))
        print('Try modifying the initial guess IC ...')
        print('  ...or modifying the number of iterations and/or the tolerance')

    if Data['flags'][2]:
        Data['IC'] = np.array([x0[0], z0, vy0[0]])
        Data['tf'] = t[-1]*2

        from .Halo_Plot import Halo_Plot
        (states_po, times_po, T_po, eigvec) = Halo_Plot(Data)

        return (states_po, times_po, T_po, eigvec)
    else:
        import os
        text = '# Data Produced by Halo Numerical Computation #\n' + '# opt = ' +\
            str(Data['opt']) + '; LP = ' + str(Data['LP']) + '; m = ' +\
            str(Data['m']) + '; phi = ' + str(Data['phi']) + '; Az = ' +\
            str(Data['Az']) + ';\n' + 'x0  = %.20f\n' % x0 + 'z0  = %.20f\n' % z0 +\
            'vy0 = %.20f\n' % vy0
        fid = open('Halo_Orbits' + os.sep + 'sample.txt','w')
        fid.write(text)
        fid.close()
