def Halo_Generator(Data):
    ##### HALO ORBITS GENERATION #####
    #
    # Importing required functions
    #
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    print('--- Halo Generator: 3rd order Richardson expansion ---')
    if Data['mode'] == 'SE':
        print('\nCR3BP: SUN-(EARTH+MOON) SYSTEM')
    else:
        print('\nCR3BP: EARTH-MOON SYSTEM')

    t  = tspan = np.linspace(0, Data['tf'], int(Data['tf']/Data['prnt_out_dt']) +1)
        # s <-- Propagation time

    Tconversion = Data['TSE']/(2*np.pi)

    # Basic Variables Calculation

    rh = (Data['mu']/3)**(1/3)
    gamma1ch = rh*(1-1/3*rh-1/9*rh**2-1/27*rh**3)
    gamma2ch = rh*(1+1/3*rh-1/9*rh**2)

    # gamma polynomial
    poly  = [1, -(3-Data['mu']), 3-2*Data['mu'], -Data['mu'], 2*Data['mu'], -Data['mu']] # for gamma 1 --> L1
    poly2 = [1, (3-Data['mu']), 3-2*Data['mu'], -Data['mu'], -2*Data['mu'], -Data['mu']] # for gamma 2 --> L2

    root1  = np.roots(poly)
    root2  = np.roots(poly2)
    gamma1 = np.real(root1[np.imag(root1) == 0])
    gamma2 = np.real(root2[np.imag(root2) == 0])

    # cn coefficients calculation
    cn1 = lambda n: 1/gamma1**3*(+1**n*Data['mu']+(-1)**n*((1-Data['mu'])*gamma1**(n+1))/(1-gamma1)**(n+1)) # L1
    cn2 = lambda n: 1/gamma2**3*((-1)**n*Data['mu']+(-1)**n*((1-Data['mu'])*gamma2**(n+1))/(1+gamma2)**(n+1)) # L2

    if Data['LP'] == 1:
        c2 = cn1(2); c3 = cn1(3); c4 = cn1(4)
        xL = (1-Data['mu'])-gamma1 # Position of Lpoint w.r.t m1-m2 center of mass
        gammaC = gamma1
    elif Data['LP'] == 2:
        c2 = cn2(2); c3 = cn2(3); c4 = cn2(4)
        xL = (1-Data['mu'])+gamma2
        gammaC = gamma2
    else:
        raise Exception('Halo Orbits:LPError.\
            The specified value for the LP parameter is outside range: [1, 2]!')

    r1 = Data['Ltheory']*((1-Data['mu'])-gamma1)
    r2 = Data['Ltheory']*((1-Data['mu'])+gamma2)

    # Eigenvalues Calculation

    wp = np.sqrt((2-c2+np.sqrt(9*c2**2-8*c2))/2)
    wv = np.sqrt(c2)
    lmbda = wp # For halo orbits

    # Variables - Linear solution

    kappa   = (wp**2+1+2*c2)/(2*wp)
    kappa_C = 2*lmbda/(lmbda**2+1-c2)

    # Third-Order Richardson Expansion

    d1 = 3*lmbda**2/kappa*(kappa*(6*lmbda**2-1)-2*lmbda)
    d2 = 8*lmbda**2/kappa*(kappa*(11*lmbda**2-1)-2*lmbda)

    a21 = (3*c3*(kappa**2-2))/(4*(1+2*c2))
    a22 = 3*c3/(4*(1+2*c2))
    a23 = -3*c3*lmbda/(4*kappa*d1)*(3*kappa**3*lmbda-6*kappa*(kappa-lmbda)+4)
    a24 = -3*c3*lmbda/(4*kappa*d1)*(2+3*kappa*lmbda)

    b21 = -3*c3*lmbda/(2*d1)*(3*kappa*lmbda-4)
    b22 = 3*c3*lmbda/d1

    d21 = -c3/(2*lmbda**2)

    a31 = -9*lmbda/(4*d2)*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa**2))+\
        (9*lmbda**2+1-c2)/(2*d2)*(3*c3*(2*a23-kappa*b21)+c4*(2+3*kappa**2))
    a32 = -1/d2*(9*lmbda/4*(4*c3*(kappa*a24-b22)+kappa*c4)+\
        3/2*(9*lmbda**2+1-c2)*(c3*(kappa*b22+d21-2*a24)-c4))

    b31 = 3/(8*d2)*(8*lmbda*(3*c3*(kappa*b21-2*a23)-c4*(2+3*kappa**2))+\
        (9*lmbda**2+1+2*c2)*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa**2)))
    b32 = 1/d2*(9*lmbda*(c3*(kappa*b22+d21-2*a24)-c4)+\
        3/8*(9*lmbda**2+1+2*c2)*(4*c3*(kappa*a24-b22)+kappa*c4))

    d31 = 3/(64*lmbda**2)*(4*c3*a24+c4)
    d32 = 3/(64*lmbda**2)*(4*c3*(a23-d21)+c4*(4+kappa**2))

    s1 = (2*lmbda*(lmbda*(1+kappa**2)-2*kappa))**-1*\
        (3/2*c3*(2*a21*(kappa**2-2)-a23*(kappa**2+2)-2*kappa*b21)-3/8*c4*(3*kappa**4-8*kappa**2+8))
    s2 = (2*lmbda*(lmbda*(1+kappa**2)-2*kappa))**-1*\
        (3/2*c3*(2*a22*(kappa**2-2)+a24*(kappa**2+2)+2*kappa*b22+5*d21)+3/8*c4*(12-kappa**2))

    l1    = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-kappa**2)+2*lmbda**2*s1
    l2    = 3/2*c3*(a24-2*a22)+9/8*c4+2*lmbda**2*s2
    incre = wp**2-wv**2

    # Relationships

    # Amplitude Constraint
    Az = Data['Az']/Data['L']
    Ax = np.sqrt(-(incre+l2*Az**2)/l1) # km

    # Phase Angle Relationship
    if Data['m'] in [1, 3]:
        psi = Data['phi']+Data['m']*np.pi/2 # rad
    else:
        raise Exception('Halo Orbits:mError.\
            The specified value for the m parameter is outside range: [1 or 3]!')

    # Added variables - Lindstedt-Poncairé Method

    nu1 = 0
    nu2 = s1*(Ax)**2+s2*(Az)**2
    nu  = 1+nu1+nu2
    tau = nu*t

    tau1   = wp*tau+Data['phi']
    deltam = 2-Data['m']

    T      = 2*np.pi/(wp*nu)*Tconversion # Period Calculation (s)
    Tad    = 2*np.pi/(wp*nu)             # Period Calculation (-)
    Tdays  = T/3600/24                   # s to days

    # EoM - 3rd Order Richardson Expansion

    xad = a21*Ax**2+a22*Az**2-Ax*np.cos(tau1)+\
        (a23*Ax**2-a24*Az**2)*np.cos(2*tau1)+(a31*Ax**3-a32*Ax*Az**2)*np.cos(3*tau1)
    yad = kappa*Ax*np.sin(tau1)+\
        (b21*Ax**2-b22*Az**2)*np.sin(2*tau1)+(b31*Ax**3-b32*Ax*Az**2)*np.sin(3*tau1)
    zad = deltam*Az*np.cos(tau1)+\
        deltam*d21*Ax*Az*(np.cos(2*tau1)-3)+deltam*(d32*Az*Ax**2-d31*Az**3)*np.cos(3*tau1)

    # Introduce rx = f(xad)   ry = f(yad)   rz = f(zad)

    x = xad*gammaC+xL; y = yad*gammaC; z = zad*gammaC # REVISAR

    # IC calculation --> For Initial Guess

    x0  = x[0];          y0 = y[0];          z0 = z[0]
    dt  = t[1]-t[0]
    vx0 = (x[1]-x0)/dt; vy0 = (y[1]-y0)/dt; vz0 = (z[1]-z0)/dt

    # # Plots

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(x, y, z)
    ax.plot(xL, 0, 0, 'ro')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    # 1 orbit display

    tdim = t*Tconversion
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, constrained_layout = True)
    ax1.plot(x[tdim<T], y[tdim<T])
    ax1.plot(0,0,'bx')
    ax1.plot(xL,0,'ro')
    ax1.set(xlabel='x', ylabel='y')
    ax2.plot(x[tdim<T], z[tdim<T])
    ax2.plot(0,0,'bx')
    ax2.plot(xL,0,'ro')
    ax2.set(xlabel='x', ylabel='z')
    ax3.plot(y[tdim<T], z[tdim<T])
    ax3.plot(0,0,'bx')
    ax3.plot(0,0,'ro')
    ax3.set(xlabel='y', ylabel='z')
    if Data['m'] == 1:
        fig.suptitle('Class I  Orbit  (L' + str(Data['LP']) +')')
    else:
        fig.suptitle('Class II  Orbit  (L' + str(Data['LP']) +')')
    plt.show()

    # # Displays
    print('-- Results ---')
    print('Ax = %.1f' % (Ax*Data['L']) + ' km;   Az = ' + str(Data['Az']) + ' km')
    print('T (days) = %.5f' % Tdays)
    if Data['m'] == 1:
        print('Case: Northern (L' + str(Data['LP']) + ')')
    else:
        print('Case: Southern (L' + str(Data['LP']) + ')')

    # IC guess display
    print('\n--- Initial Guess Display ---')
    print('x0  = %.20f;' % x0)
    print('y0  = %.20f;' % y0)
    print('z0  = %.20f;' % z0)
    print('vx0 = %.20f;' % vx0)
    print('vy0 = %.20f;' % vy0)
    print('vz0 = %.20f;\n' % vz0)

    # Time guesses
    print('--- Time Guess ---')
    print('T  = %.20f;' % Tad)
    print('T/2  = %.20f;\n' % (Tad/2))

    if Data['flags'][1]:
        Data['IC'] = np.array([x0, z0, vy0])

        from .Halo_Num_Comp import Halo_Num_Comp
        (states_po, T_po, eigvec) = Halo_Num_Comp(Data)

        return (states_po, T_po, eigvec)
    else:
        import os
        text = '# Data Produced by Halo Generator #\n' + '# mode = ' + str(Data['mode']) +\
            '; LP = ' + str(Data['LP']) + '; m = ' + str(Data['m']) + '; phi = ' +\
            str(Data['phi']) + '; Az = ' + str(Data['Az']) + ';\n' +\
            'x0  = %.20f\n' % x0 + 'y0  = %.20f\n' % y0 + 'z0  = %.20f\n' % z0 +\
            'vx0 = %.20f\n' % vx0 + 'vy0 = %.20f\n' % vy0 + 'vz0 = %.20f\n' % vz0
        fid = open('Halo_Orbits' + os.sep + 'sample.txt','w')
        fid.write(text)
        fid.close()
