def construct(params, T0, states_po, times_po, eigvec, eigval, inv_phi_0,
    prnt_out_dt, npoints, stop_fun, ang, L):

    #
    # Import required functions
    #
    import numpy as np
    from scipy import linalg
    from PCR3BP import PCR3BP_propagator

    ## Step 4: Computation of invariant manifolds
    T_po = T0
    idx  = np.linspace(0, len(times_po), npoints, endpoint = False, dtype = int)
    et0  = 0
    eps  = 1e-7
    deltat = 5*T_po
    sign   = np.array([-1, 1])

    # Matrix initialization
    states_u = []
    SF_u     = np.array([])
    times_u  = []
    states_s = []
    SF_s     = np.array([])
    times_s  = []

    for i in range(npoints):
        for j in range(len(sign)):
            x = states_po[:, idx[i]]
            t = times_po[idx[i]]
            phi = [eigvec[:, 0]*np.exp(eigval[0]*t), eigvec[:, 1]*np.exp(eigval[1]*t),
                eigvec[:, 2]*np.exp(eigval[2]*t), eigvec[:, 3]*np.exp(eigval[3]*t)]
            Phi = phi * inv_phi_0
            [mon_eigvals, mon_eigvecs] = linalg.eig(Phi)

            Yu = mon_eigvecs[:, 0] # unstable eigenvector
            Ys = mon_eigvecs[:, 1] # stable eigenvector

            Xu = np.real(states_po[:, idx[i]] + eps*sign[j]*Yu)
            Xs = np.real(states_po[:, idx[i]] + eps*sign[j]*Ys)

            # Integrate unstable manifold forwards in time
            [SF, etf_u, states, times] = PCR3BP_propagator(Xu, et0, deltat,
                prnt_out_dt, stop_fun, params[0], params[1], ang, L)

            states_u.append(states)
            times_u.append(times)
            SF_u = np.append([SF_u], [SF])

            # Integrate stable manifold backwards in time
            [SF, etf, states, times] =  PCR3BP_propagator(Xs, et0, -deltat,
                prnt_out_dt, stop_fun, params[0], params[1], ang, L)

            states_s.append(states)
            times_s.append(times)
            SF_s = np.append([SF_s], [SF])

            print('Iteration: ' + str(2*i + j +1))

    return [states_s, times_s, SF_s.reshape(-1, 4), states_u, times_u, SF_u.reshape(-1, 4)]

def plotm(mu1, mu2, pos, states_po, states_s, SF_s, states_u, SF_u, ang):

    import matplotlib.pyplot as plt
    import numpy as np

    plt.plot(mu1, 0, 'ko')
    plt.plot(pos[1, 0], pos[1, 1], 'ko')
    plt.plot(states_po[0], states_po[1], 'k')
    rmax = (states_po[0] - pos[1, 0])**2 + (states_po[1] - pos[1, 1])**2\
        == max((states_po[0] - pos[1, 0])**2 + (states_po[1] - pos[1, 1])**2)
    tol  = 5e-3
    plt.plot([mu1, pos[1, 0]], [0, pos[1,0]*(states_po[1][rmax][0] - tol)\
        /(states_po[0][rmax][0] - tol)], 'g--')
    plt.plot([mu1, pos[1, 0]], [0, -pos[1,0]*(states_po[1][rmax][0] - tol)\
        /(states_po[0][rmax][0] - tol)], 'g--')
    plt.plot([mu1, mu1 + (pos[1, 0] - mu1)*np.cos(ang*np.pi/180)],
        [0, (pos[1, 0] - mu1)*np.sin(ang*np.pi/180)], 'r--')
    plt.show()

    for i in range(len(states_u)):
        plt.plot(states_u[i][0], states_u[i][1], 'r')
        plt.plot(states_s[i][0], states_s[i][1], 'b')
    plt.plot(mu1, 0 , 'ro')
    plt.plot(pos[1, 0], pos[1, 1], 'ko')
    plt.plot(states_po[0], states_po[1], 'k')
    plt.plot([mu1, mu1 + (pos[1, 0] - mu1)*np.cos(ang*np.pi/180)],
        [0, (pos[1, 0] - mu1)*np.sin(ang*np.pi/180)], 'y--')
    plt.show()

    plt.plot(mu1, 0 , 'ro')
    plt.plot(pos[:2, 0], pos[:2, 1], 'ko')
    plt.show()

    tangVels = max(abs(np.cos(ang*np.pi/180)*SF_s[SF_s[:, 1] < (pos[1, 0]-mu1), 2]\
        - np.sin(ang*np.pi/180)*SF_s[SF_s[:, 1] < (pos[1, 0]-mu1), 3]))
    tangVelu = max(abs(np.cos(ang*np.pi/180)*SF_u[SF_u[:, 1] > -(pos[1, 0]-mu1), 2]\
        - np.sin(ang*np.pi/180)*SF_u[SF_u[:, 1] > -(pos[1, 0]-mu1), 3]))
    print('The maximum tangential velocities registered in this Poincaré section are: ')
    print('Max Vt (stable manifold):   %1.5f' % tangVels)
    print('Max Vt (unstable manifold): %1.5f' % tangVelu)

    fig, ax = plt.subplots()
    for j in range(len(SF_s)):
        if SF_s[j, 1] < (pos[1, 0]-mu1):
            ax.plot(np.sqrt((SF_s[j, 0] - mu1)**2 + SF_s[j, 1]**2),
                np.sin(ang*np.pi/180)*SF_s[j, 2] + np.cos(ang*np.pi/180)*SF_s[j, 3],
                'bo')
        if SF_u[j, 1] > -(pos[1, 0]-mu1):
            ax.plot(np.sqrt((SF_u[j, 0] - mu1)**2 + SF_u[j, 1]**2),
                np.sin(ang*np.pi/180)*SF_u[j, 2] + np.cos(ang*np.pi/180)*SF_u[j, 3],
                'ro')
    xu = np.sqrt((SF_u[SF_u[:, 1] > -(pos[1, 0]-mu1), 0] - mu1)**2\
        + SF_u[SF_u[:, 1] > -(pos[1, 0]-mu1), 1]**2)
    yu = np.sin(ang*np.pi/180)*SF_u[SF_u[:, 1] > -(pos[1, 0]-mu1), 2]\
        + np.cos(ang*np.pi/180)*SF_u[SF_u[:, 1] > -(pos[1, 0]-mu1), 3]
    xs = np.sqrt((SF_s[SF_s[:, 1] < (pos[1, 0]-mu1), 0] - mu1)**2\
        + SF_s[SF_s[:, 1] < (pos[1, 0]-mu1), 1]**2)
    ys = np.sin(ang*np.pi/180)*SF_s[SF_s[:, 1] < (pos[1, 0]-mu1), 2]\
        + np.cos(ang*np.pi/180)*SF_s[SF_s[:, 1] < (pos[1, 0]-mu1), 3]

    centeru = [np.mean(xu), np.mean(yu)]
    centers = [np.mean(xs), np.mean(ys)]
    atanu   = np.arctan2(yu - centeru[1], xu - centeru[0])
    atans   = np.arctan2(ys - centers[1], xs - centers[0])
    keyu    = np.argsort(atanu)
    keys    = np.argsort(atans)

    ax.fill(xu[keyu], yu[keyu], facecolor = 'salmon')
    ax.fill(xs[keys], ys[keys], facecolor = 'lightblue')
    plt.xlabel('Distance from 2nd primary')
    plt.ylabel('Velocity normal to the plane')
    plt.title(r'Poincaré section at angle %2.1f$^{\circ}$' % ang)
    plt.show()
