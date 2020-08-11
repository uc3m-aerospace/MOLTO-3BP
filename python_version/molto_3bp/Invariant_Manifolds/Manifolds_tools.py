def construct(params, T0, states_po, times_po, eigvec, eigval, inv_phi_0,
    prnt_out_dt, npoints, stop_fun, ang, L):

    #
    # Import required functions
    #
    import numpy as np
    from scipy import linalg
    from PCR3BP import PCR3BP_propagator

    ## Computation of invariant manifolds
    T_po = T0
    idx  = np.linspace(0, len(times_po), npoints, endpoint = False, dtype = int)
    et0  = 0
    eps  = 1e-7
    deltat = 5*T_po
    sign   = np.array([-1, 1])
    veclen = len(states_po[:, 0])

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
            t = times_po[idx[i]]/100
            phi = np.zeros((veclen, veclen)) + 0.j
            for k in range(veclen):
                phi[:, k] = eigvec[:, k]*np.exp(eigval[k]*t)
            Phi = phi * inv_phi_0
            [mon_eigvals, mon_eigvecs] = linalg.eig(Phi)

            Yu = mon_eigvecs[:, 0] # unstable eigenvector
            Ys = mon_eigvecs[:, 1] # stable eigenvector

            Xu = np.real(states_po[:, idx[i]] + eps*sign[j]*Yu)
            Xs = np.real(states_po[:, idx[i]] + eps*sign[j]*Ys)

            stop_fun.direction = -(L - params[0])

            # Integrate unstable manifold forwards in time
            [SF, etf_u, states, times] = PCR3BP_propagator(Xu, et0, deltat,
                prnt_out_dt, stop_fun, params[0], params[1], ang, L)

            states_u.append(states)
            times_u.append(times)
            SF_u = np.append([SF_u], [SF])

            stop_fun.direction = (L - params[0])

            # Integrate stable manifold backwards in time
            [SF, etf, states, times] =  PCR3BP_propagator(Xs, et0, -deltat,
                prnt_out_dt, stop_fun, params[0], params[1], ang, L)

            states_s.append(states)
            times_s.append(times)
            SF_s = np.append([SF_s], [SF])

            print('Iteration: ' + str(2*i + j +1))

    return [states_s, times_s, SF_s.reshape(-1, veclen), states_u, times_u,
        SF_u.reshape(-1, veclen)]

def plotm(mu1, mu2, pos, states_po, states_s, SF_s, states_u, SF_u, ang, angmin):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    plt.plot(mu1, 0, 'ko')
    plt.plot(pos[0], pos[1], 'ko')
    plt.plot(states_po[0], states_po[1], 'k')

    plt.plot([mu1, mu1 + 1.25*abs(pos[0] - mu1)*np.cos(angmin)],
        [0, 1.25*abs(pos[0] - mu1)*np.sin(angmin)], 'g--')
    plt.plot([mu1, mu1 + 1.25*abs(pos[0] - mu1)*np.cos(angmin)],
        [0, -1.25*abs(pos[0] - mu1)*np.sin(angmin)], 'g--')
    plt.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
        [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)], 'r--')
    plt.show()

    if len(SF_s[0]) == 4:

        for i in range(len(states_u)):
            plt.plot(states_u[i][0], states_u[i][1], 'r')
            plt.plot(states_s[i][0], states_s[i][1], 'b')
        plt.plot(mu1, 0 , 'ro')
        plt.plot(pos[0], pos[1], 'ko')
        plt.plot(states_po[0], states_po[1], 'k')
        plt.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)], 'y--')
        plt.show()

        tangVels = max(abs(np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 2]\
            + np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]))
        tangVelu = max(abs(np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 2]\
            + np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]))
        print('The maximum tangential velocities registered in this Poincaré section are: ')
        print('Max Vt (stable manifold):   %1.5f' % tangVels)
        print('Max Vt (unstable manifold): %1.5f' % tangVelu)

        fig, ax = plt.subplots()
        for j in range(len(SF_s)):
            if abs(SF_s[j, 1]) < abs(pos[0]-mu1):
                ax.plot(np.sqrt((SF_s[j, 0] - mu1)**2 + SF_s[j, 1]**2),
                    - np.sin(ang*np.pi/180)*SF_s[j, 2] + np.cos(ang*np.pi/180)*SF_s[j, 3],
                    'bo')
            if abs(SF_u[j, 1]) < abs(pos[0]-mu1):
                ax.plot(np.sqrt((SF_u[j, 0] - mu1)**2 + SF_u[j, 1]**2),
                    - np.sin(ang*np.pi/180)*SF_u[j, 2] + np.cos(ang*np.pi/180)*SF_u[j, 3],
                    'ro')
        xu = np.sqrt((SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
            + SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 1]**2)
        yu = - np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 2]\
            + np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]
        xs = np.sqrt((SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
            + SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 1]**2)
        ys = - np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 2]\
            + np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]

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

    else:

        fig = plt.figure()
        ax = Axes3D(fig)
        for i in range(len(states_u)):
            ax.plot(states_u[i][0], states_u[i][1], states_u[i][2], 'r')
            ax.plot(states_s[i][0], states_s[i][1], states_s[i][2], 'b')
        ax.plot([mu1], [0], [0], 'ro')
        ax.plot([pos[0]], [pos[1]], [0], 'ko')
        ax.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)],
            [abs(pos[0] - mu1)/2, abs(pos[0] - mu1)/2], 'y--')
        ax.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)],
            [-abs(pos[0] - mu1)/2, -abs(pos[0] - mu1)/2], 'y--')
        ax.plot([mu1, mu1], [0, 0],
            [-abs(pos[0] - mu1)/2, abs(pos[0] - mu1)/2], 'y--')
        ax.plot([mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180),
            mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [abs(pos[0] - mu1)*np.sin(ang*np.pi/180),
            abs(pos[0] - mu1)*np.sin(ang*np.pi/180)],
            [-abs(pos[0] - mu1)/2, abs(pos[0] - mu1)/2], 'y--')
        ax.set(xlabel = 'x', ylabel = 'y', zlabel = 'z')
        plt.show()

        tangVels1 = max(abs(np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]\
            + np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 4]))
        tangVels2 = max(abs(SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 5]))

        tangVelu1 = max(abs(np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]\
            + np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 4]))
        tangVelu2 = max(abs(SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 5]))

        print('The maximum tangential velocities registered in this Poincaré section are: ')
        print('Max Vt1 (stable manifold):   %1.5f' % tangVels1)
        print('Max Vt2 (stable manifold):   %1.5f' % tangVels2)
        print('Max Vt1 (unstable manifold): %1.5f' % tangVelu1)
        print('Max Vt2 (unstable manifold): %1.5f' % tangVelu2)

        fig, ax = plt.subplots()
        for j in range(len(SF_s)):
            if abs(SF_s[j, 1]) < abs(pos[0]-mu1):
                ax.plot(np.sqrt((SF_s[j, 0] - mu1)**2 + SF_s[j, 1]**2),
                    - np.sin(ang*np.pi/180)*SF_s[j, 3] + np.cos(ang*np.pi/180)*SF_s[j, 4],
                    'bo')
            if abs(SF_u[j, 1]) < abs(pos[0]-mu1):
                ax.plot(np.sqrt((SF_u[j, 0] - mu1)**2 + SF_u[j, 1]**2),
                    - np.sin(ang*np.pi/180)*SF_u[j, 3] + np.cos(ang*np.pi/180)*SF_u[j, 4],
                    'ro')
        xu = np.sqrt((SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
            + SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 1]**2)
        yu = - np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]\
            + np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 4]
        xs = np.sqrt((SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
            + SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 1]**2)
        ys = - np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]\
            + np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 4]

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
