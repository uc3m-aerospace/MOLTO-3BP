def construct(params, Orbit, prnt_out_dt, npoints, d, branch, stop_fun, ang, L):

    #
    # Import required functions
    #
    import numpy as np
    from scipy import linalg
    from PCR3BP import PCR3BP_propagator

    ## Computation of invariant manifolds
    idx  = np.linspace(0, len(Orbit[0][0]), npoints, endpoint = False, dtype = int)
    et0  = 0
    deltat = 5 * Orbit[1]
    veclen = len(Orbit[0][:, 0])

    if not d:
        sign = np.array([-1, 1])
    elif d == 1 or d == -1:
        if veclen == 6 and L > params[0]:
            sign = np.array([d])
        else:
            sign = np.array([-d])
    else:
        raise Exception('Manifolds_tools:dError.'+\
            '    The direction selected is not valid [1, 0, -1]!')

    # Matrix initialization
    states_u = []
    SF_u     = np.array([])
    times_u  = []
    states_s = []
    SF_s     = np.array([])
    times_s  = []

    if not branch:
        u = 1
        s = 1
    else:
        if branch == 1:
            s = 1
            u = 0
        elif branch == -1:
            s = 0
            u = 1
        else:
            raise Exception('Manifolds_tools:branchError.'+\
                '    The branch selected is not valid [1, 0, -1]!')

    if veclen == 4:
        eps = 1e-7
    else:
        eps = 1e-4

    for i in range(npoints):
        for j in range(len(sign)):

            print('Iteration: ' + str(len(sign)*i + j +1))

            if u:
                Yu = Orbit[2][:, 0] # unstable eigenvector

                if veclen == 6 or L > params[0]:
                    Xu = np.real(Orbit[0][:, idx[i]] + eps*sign[j]*Yu)
                else:
                    Xu = np.real(Orbit[0][:, idx[i]] - eps*sign[j]*Yu)

                stop_fun.direction = -(L - params[0])

                # Integrate unstable manifold forwards in time
                [SF, etf_u, states, times] = PCR3BP_propagator(Xu, et0, deltat,
                    prnt_out_dt, stop_fun, params[0], params[1], ang, L)

                states_u.append(states)
                times_u.append(times)
                SF_u = np.append([SF_u], [SF])

            if s:
                Ys = Orbit[2][:, 1] # stable eigenvector

                Xs = np.real(Orbit[0][:, idx[i]] - eps*sign[j]*Ys)

                stop_fun.direction = (L - params[0])

                # Integrate stable manifold backwards in time
                [SF, etf, states, times] =  PCR3BP_propagator(Xs, et0, -deltat,
                    prnt_out_dt, stop_fun, params[0], params[1], ang, L)

                states_s.append(states)
                times_s.append(times)
                SF_s = np.append([SF_s], [SF])

    return [states_s, times_s, SF_s.reshape(-1, veclen), states_u, times_u,
        SF_u.reshape(-1, veclen)]

def fourierTest(mu1, mu2, pos, states_s, states_u, ang, Data):

    import numpy as np
    from scipy.fft import fft, fftfreq, fftshift
    import matplotlib.pyplot as plt

    signal_s = []
    signal_u = []

    fig1, ax1   = plt.subplots()
    fig12, ax12 = plt.subplots()
    fig2, ax2   = plt.subplots()
    fig22, ax22   = plt.subplots()
    fig3, ax3   = plt.subplots()
    fig32, ax32   = plt.subplots()

    if len(states_s):
        for i in states_s:
            if abs(i[-1, 1]) < abs(pos[0] - mu1):
                if len(i[0, :0:-1]) % 2:
                    temp = fft(i[0, :1:-1]-pos[0] + i[1, :1:-1]*1.j)
                    l = len(i[0, :1:-1])
                    signal_s.append(temp)
                else:
                    temp = fft(i[0, :0:-1]-pos[0] + i[1, :0:-1]*1.j)
                    l = len(i[0, :0:-1])
                    signal_s.append(temp)
                xf = fftfreq(l, Data['prnt_out_dt'])
                xf = fftshift(xf)
                yf = fftshift(temp)
                ax1.plot(xf[abs(yf) > 0.01*max(abs(yf))],
                    100.0/l*np.real(yf[abs(yf) > 0.01*max(abs(yf))]))
                ax12.plot(xf[abs(yf) > 0.01*max(abs(yf))],
                    100.0/l*np.imag(yf[abs(yf) > 0.01*max(abs(yf))]))

    if len(states_u):
        for i in states_u:
            if abs(i[-1, 1]) < abs(pos[0] - mu1):
                if len(i[0, :-1]) % 2:
                    temp = fft(i[0, :-2]-pos[0] + i[1, :-2]*1.j)
                    l = len(i[0, :-2])
                    signal_u.append(temp)
                else:
                    temp = fft(i[0, :-1]-pos[0] + i[1, :-1]*1.j)
                    l = len(i[0, :-1])
                    signal_u.append(temp)
                xf = fftfreq(l, Data['prnt_out_dt'])
                xf = fftshift(xf)
                yf = fftshift(temp)
                ax2.plot(xf[abs(yf) > 0.01*max(abs(yf))],
                    100.0/l*np.real(yf[abs(yf) > 0.01*max(abs(yf))]))
                ax22.plot(xf[abs(yf) > 0.01*max(abs(yf))],
                    100.0/l*np.imag(yf[abs(yf) > 0.01*max(abs(yf))]))

    if len(states_u) and len(states_s) and Data['d'] == 1:
        for i in range(len(states_s)):
            if len(states_s[i][0, :-1]) % 2:
                vecx1 = states_s[i][0, :1:-1]
                vecy1 = states_s[i][1, :1:-1]
            else:
                vecx1 = states_s[i][0, :0:-1]
                vecy1 = states_s[i][1, :0:-1]

            if len(states_u[i][0, :-1]) % 2:
                vecx2 = states_u[i][0, :-2]
                vecy2 = states_u[i][1, :-2]
            else:
                vecx2 = states_u[i][0, :-1]
                vecy2 = states_u[i][1, :-1]

            vecx = np.append(vecx1, vecx2)
            vecy = np.append(vecy1, vecy2)

            temp = fft(vecx-pos[0] + vecy*1.j)
            l = len(vecx)

            signal_s.append(temp)

            xf = fftfreq(l, Data['prnt_out_dt'])
            xf = fftshift(xf)
            yf = fftshift(temp)
            ax3.plot(xf[abs(yf) > 0.01*max(abs(yf))],
                100.0/l*np.real(yf[abs(yf) > 0.01*max(abs(yf))]))
            ax32.plot(xf[abs(yf) > 0.01*max(abs(yf))],
                100.0/l*np.imag(yf[abs(yf) > 0.01*max(abs(yf))]))

    plt.show()

def plotm(mu1, mu2, pos, states_po, states_s, SF_s, states_u, SF_u, ang, angmin):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    fig1, ax1 = plt.subplots()
    ax1.plot(mu1, 0, 'ko')
    ax1.plot(pos[0], pos[1], 'ko')

    if len(states_po) == 2:
        ax1.plot(states_po[0][0][0], states_po[0][0][1], 'k')
        ax1.plot(states_po[1][0][0], states_po[1][0][1], 'k')
        ax1.plot(pos[2], pos[3], 'ko')
        size = len(states_po[0][0][:, 0])
    else:
        ax1.plot(states_po[0][0], states_po[0][1], 'k')
        size = len(states_po[0][:, 0])

    ax1.plot([mu1, mu1 + 1.25*abs(pos[0] - mu1)*np.cos(angmin)],
        [0, 1.25*abs(pos[0] - mu1)*np.sin(angmin)], 'g--')
    ax1.plot([mu1, mu1 + 1.25*abs(pos[0] - mu1)*np.cos(angmin)],
        [0, -1.25*abs(pos[0] - mu1)*np.sin(angmin)], 'g--')
    ax1.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
        [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)], 'r--')

    if size == 4:

        fig2, ax2 = plt.subplots()
        for i_u in range(len(states_u)):
            ax2.plot(states_u[i_u][0], states_u[i_u][1], 'r')
        for i_s in range(len(states_s)):
            ax2.plot(states_s[i_s][0], states_s[i_s][1], 'b')
        ax2.plot(mu1, 0 , 'ro')
        ax2.plot(pos[0], pos[1], 'ko')
        if len(states_po) == 2:
            ax1.plot(states_po[0][0][0], states_po[0][0][1], 'k')
            ax1.plot(states_po[1][0][0], states_po[1][0][1], 'k')
            ax1.plot(pos[2], pos[3], 'ko')
            size = len(states_po[0][0][:, 0])
        else:
            ax1.plot(states_po[0][0], states_po[0][1], 'k')
            size = len(states_po[0][:, 0])
        ax2.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)], 'y--')

        print('The maximum tangential velocities registered in this Poincaré section are: ')

        u = 0
        s = 0
        if len(SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), :]) > 0:
            tangVels = max(abs(np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 2]\
                + np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]))
            print('Max Vt (stable manifold):   %1.5f' % tangVels)

            xs = np.sqrt((SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
                + SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 1]**2)
            ys = - np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 2]\
                + np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]

            centers = [np.mean(xs), np.mean(ys)]
            atans   = np.arctan2(ys - centers[1], xs - centers[0])
            keys    = np.argsort(atans)

            s = 1

        if len(SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), :]) > 0:
            tangVelu = max(abs(np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 2]\
                + np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]))
            print('Max Vt (unstable manifold): %1.5f' % tangVelu)

            xu = np.sqrt((SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
                + SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 1]**2)
            yu = - np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 2]\
                + np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]

            centeru = [np.mean(xu), np.mean(yu)]
            atanu   = np.arctan2(yu - centeru[1], xu - centeru[0])
            keyu    = np.argsort(atanu)

            u = 1

        if u or s:
            fig3, ax3 = plt.subplots()
            for j in range(len(SF_s)):
                if abs(SF_s[j, 1]) < abs(pos[0]-mu1):
                    ax3.plot(np.sqrt((SF_s[j, 0] - mu1)**2 + SF_s[j, 1]**2),
                        - np.sin(ang*np.pi/180)*SF_s[j, 2] + np.cos(ang*np.pi/180)*SF_s[j, 3],
                        'bo')
            for k in range(len(SF_u)):
                if abs(SF_u[k, 1]) < abs(pos[0]-mu1):
                    ax3.plot(np.sqrt((SF_u[k, 0] - mu1)**2 + SF_u[k, 1]**2),
                        - np.sin(ang*np.pi/180)*SF_u[k, 2] + np.cos(ang*np.pi/180)*SF_u[k, 3],
                        'ro')

            ax3.set(xlabel = 'Distance from 2nd primary')
            ax3.set(ylabel = 'Velocity normal to the plane')
            ax3.set(title = r'Poincaré section at angle %2.1f$^{\circ}$' % ang)

        if u:
            ax3.fill(xu[keyu], yu[keyu], facecolor = 'salmon')
        if s:
            ax3.fill(xs[keys], ys[keys], facecolor = 'lightblue')

        plt.show()

    else:

        fig2 = plt.figure()
        ax2 = Axes3D(fig2)
        for i_u in range(len(states_u)):
            ax2.plot(states_u[i_u][0], states_u[i_u][1], states_u[i_u][2], 'r')
        for i_s in range(len(states_s)):
            ax2.plot(states_s[i_s][0], states_s[i_s][1], states_s[i_s][2], 'b')
        ax2.plot([mu1], [0], [0], 'ro')
        ax2.plot([pos[0]], [pos[1]], [0], 'ko')
        ax2.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)],
            [abs(pos[0] - mu1)/2, abs(pos[0] - mu1)/2], 'y--')
        ax2.plot([mu1, mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [0, abs(pos[0] - mu1)*np.sin(ang*np.pi/180)],
            [-abs(pos[0] - mu1)/2, -abs(pos[0] - mu1)/2], 'y--')
        ax2.plot([mu1, mu1], [0, 0],
            [-abs(pos[0] - mu1)/2, abs(pos[0] - mu1)/2], 'y--')
        ax2.plot([mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180),
            mu1 + abs(pos[0] - mu1)*np.cos(ang*np.pi/180)],
            [abs(pos[0] - mu1)*np.sin(ang*np.pi/180),
            abs(pos[0] - mu1)*np.sin(ang*np.pi/180)],
            [-abs(pos[0] - mu1)/2, abs(pos[0] - mu1)/2], 'y--')
        ax2.set(xlabel = 'x', ylabel = 'y', zlabel = 'z')

        print('The maximum tangential velocities registered in this Poincaré section are: ')

        u = 0
        s = 0

        if len(SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), :]) > 0:
            tangVels1 = max(abs(np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]\
                + np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 4]))
            print('Max Vt1 (stable manifold):   %1.5f' % tangVels1)
            tangVels2 = max(abs(SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 5]))
            print('Max Vt2 (stable manifold):   %1.5f' % tangVels2)

            xs = np.sqrt((SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
                + SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 1]**2\
                + SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 2]**2)
            ys = - np.sin(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 3]\
                + np.cos(ang*np.pi/180)*SF_s[abs(SF_s[:, 1]) < abs(pos[0]-mu1), 4]

            centers = [np.mean(xs), np.mean(ys)]
            atans   = np.arctan2(ys - centers[1], xs - centers[0])
            keys    = np.argsort(atans)

            s = 1

        if len(SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), :]) > 0:
            tangVelu1 = max(abs(np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]\
                + np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 4]))
            print('Max Vt1 (unstable manifold): %1.5f' % tangVelu1)
            tangVelu2 = max(abs(SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 5]))
            print('Max Vt2 (unstable manifold): %1.5f' % tangVelu2)

            xu = np.sqrt((SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 0] - mu1)**2\
                + SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 1]**2\
                + SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 2]**2)
            yu = - np.sin(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 3]\
                + np.cos(ang*np.pi/180)*SF_u[abs(SF_u[:, 1]) < abs(pos[0]-mu1), 4]

            centeru = [np.mean(xu), np.mean(yu)]
            atanu   = np.arctan2(yu - centeru[1], xu - centeru[0])
            keyu    = np.argsort(atanu)

            u = 1

        if u or s:
            fig3, ax3 = plt.subplots()
            for j in range(len(SF_s)):
                if abs(SF_s[j, 1]) < abs(pos[0]-mu1):
                    ax3.plot(np.sqrt((SF_s[j, 0] - mu1)**2 + SF_s[j, 1]**2\
                        + SF_s[j, 2]**2), -np.sin(ang*np.pi/180)*SF_s[j, 3]\
                        + np.cos(ang*np.pi/180)*SF_s[j, 4],
                        'bo')
            for k in range(len(SF_u)):
                if abs(SF_u[k, 1]) < abs(pos[0]-mu1):
                    ax3.plot(np.sqrt((SF_u[k, 0] - mu1)**2 + SF_u[k, 1]**2\
                        + SF_u[k, 2]**2), -np.sin(ang*np.pi/180)*SF_u[k, 3]\
                        + np.cos(ang*np.pi/180)*SF_u[k, 4],
                        'ro')

            ax3.set(xlabel = 'Distance from 2nd primary')
            ax3.set(ylabel = 'Velocity normal to the plane')
            ax3.set(title  = r'Poincaré section at angle %2.1f$^{\circ}$' % ang)

            if s:
                ax3.fill(xs[keys], ys[keys], facecolor = 'lightblue')
            if u:
                ax3.fill(xu[keyu], yu[keyu], facecolor = 'salmon')

        plt.show()
