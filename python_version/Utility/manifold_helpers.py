import Utility.helpers as h
import Utility.exceptions as exe
import Common.Config.configurator as conf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import Utility.lagrange_helpers as lh
import numpy as np
from scipy import linalg
from Utility.PCR3BP_helpers import PCR3BP_propagator
import os


def query_func():
    h.log('Mission sequence simulation')

    # Iteration variables preallocated
    seq = {'it': 1}

    while 1:

        h.log('\nParameters for Orbit ' + str(seq['it']))

        seq['type' + str(seq['it'])] = input('Orbit type [HL, LY]: ')

        if seq['type' + str(seq['it'])] == 'LY':
            Ax = input('Amplitude (Ax, nondimensional): ')
            LP = input('Lagrange Point [1, 2]: ')

            seq['Ax' + str(seq['it'])] = float(Ax)
            seq['LP' + str(seq['it'])] = int(LP)

        elif seq['type' + str(seq['it'])] == 'HL':
            Az = input('Amplitude (Az, dimensional [km]): ')
            phi = input('Azimuthal rotation [rads]: ')
            m = input('North/South variants of the orbit [1, 3]: ')
            LP = input('Lagrange Point [1, 2]: ')

            seq['Az' + str(seq['it'])] = float(Az)
            seq['phi' + str(seq['it'])] = float(phi)
            seq['m' + str(seq['it'])] = int(m)
            seq['LP' + str(seq['it'])] = int(LP)

        else:
            config = conf.Configurator()
            if config.get('run.mode', 'debug').lower() != 'debug':
                raise exe.ManifoldsSequenceError()
            else:
                exe.ManifoldsSequenceError()

        key = input('\nWould you like to add a section? (1 == yes, 0 == no): ')

        if not int(key):
            break

        h.log('\nParameters for Section ' + str(seq['it']))
        ang = input('Angle between section and +X semiplane [degrees]: ')

        seq['ang' + str(seq['it'])] = int(ang)

        key = input('\nWould you like to add another orbit? (1 == yes, 0 == no): ')

        if not int(key):
            break

        seq['it'] = seq['it'] + 1

    return seq


def query_func_text(Input):
    h.log('Mission sequence simulation')

    # Iteration variables preallocated
    seq = {'it': 1}
    fid = open(os.getcwd() + Input['file'], 'r')
    data = fid.readlines()
    refdata = []
    for i in range(len(data)):
        if data[i].split()[0][0] != '#':
            refdata.append(i)

    data = [data[j] for j in refdata]
    line = 0

    while 1:

        seq['type' + str(seq['it'])] = data[line].split()[-1]

        if seq['type' + str(seq['it'])] == 'LY':
            Ax = data[line + 1].split()[-1]
            LP = data[line + 2].split()[-1]

            seq['Ax' + str(seq['it'])] = float(Ax)
            seq['LP' + str(seq['it'])] = int(LP)

            if len(data) == line + 3:
                break

            ang = data[line + 3].split()[-1]

            seq['ang' + str(seq['it'])] = int(ang)

            if len(data) == line + 4:
                break

            line = line + 4

        elif seq['type' + str(seq['it'])] == 'HL':
            Az = data[line + 1].split()[-1]
            phi = data[line + 2].split()[-1]
            m = data[line + 3].split()[-1]
            LP = data[line + 4].split()[-1]

            seq['Az' + str(seq['it'])] = float(Az)
            seq['phi' + str(seq['it'])] = float(phi)
            seq['m' + str(seq['it'])] = int(m)
            seq['LP' + str(seq['it'])] = int(LP)

            if len(data) == line + 5:
                break

            ang = data[line + 5].split()[-1]

            seq['ang' + str(seq['it'])] = int(ang)

            if len(data) == line + 6:
                break

            line = line + 6

        else:
            config = conf.Configurator()
            if config.get('run.mode', 'debug').lower() != 'debug':
                raise exe.ManifoldsSequenceError()
            else:
                exe.ManifoldsSequenceError()

        seq['it'] = seq['it'] + 1

    fid.close()

    return seq


### -- Will load the init variables -- ####
def load_variables(data: dict):
    #
    # Import required functions
    #
    import numpy as np
    import spiceypy as spice

    ## Load Spice kernels
    h.load_kernel()

    ## INITIAL DATA

    # Gravitational constants:
    mu_earth = spice.bodvrd('EARTH', 'GM', 1)[1][0]  # Earth grav. pot [km3/s2]
    mu_sun = spice.bodvrd('SUN', 'GM', 1)[1][0]  # Sun grav. pot [km3/s2]
    mu_moon = spice.bodvrd('MOON', 'GM', 1)[1][0]  # Moon grav. pot [km3/s2]

    if data.get('mode', 'SE') == 'SE':  # Earth-Moon system
        data['mu'] = (mu_earth + mu_moon) / (mu_earth + mu_sun + mu_moon)
        # Normalization relations
        data['L'] = 1.497610041e6  # km        # Distance Sun-Earth
        data['Ltheory'] = 1.496e8
        data['TSE'] = 3.155815e7  # s            # Period of Earth around Sun: 365 days
    elif data.get('mode', 'EM') == 'EM':  # Sun-Earth system
        data['mu'] = mu_moon / (mu_moon + mu_earth)
        data['L'] = 3.84388174e5  # distance between primaries
        data['Ltheory'] = 3.84388174e5
        data['TSE'] = 2.361e6  # s               # orbital period of primaries 27.32 days
    else:
        raise exe.HaloOrbitsmodeError()

    data['params'] = (1 - data['mu'], data['mu'])

    guess = np.array([0.9, 1.01, -1])

    # Returns, x and y coordinates of L points
    data['pos'] = lh.lagrange_points(data['params'][1], guess)

    # Plot Lagrange points
    lh.plot_lagrange_points(data['params'][0], data['params'][1], data['pos'])

    return data


### --------- Tools --------

def construct(params, Orbit, prnt_out_dt, npoints, d, branch, stop_fun, ang, L):
    #
    # Import required functions
    #

    ## Computation of invariant manifolds
    idx = np.linspace(0, len(Orbit[0][0]), npoints, endpoint=False, dtype=int)
    et0 = 0
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
        raise Exception('Manifolds_tools:dError.' + \
                        '    The direction selected is not valid [1, 0, -1]!')

    # Matrix initialization
    states_u = []
    SF_u = np.array([])
    times_u = []
    states_s = []
    SF_s = np.array([])
    times_s = []

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
            raise Exception('Manifolds_tools:branchError.' + \
                            '    The branch selected is not valid [1, 0, -1]!')

    if veclen == 4:
        eps = 1e-7
    else:
        eps = 1e-4

    for i in range(npoints):
        for j in range(len(sign)):

            h.log('Iteration: ' + str(len(sign) * i + j + 1))

            if u:
                Yu = Orbit[2][:, 0]  # unstable eigenvector

                if veclen == 6 or L > params[0]:
                    Xu = np.real(Orbit[0][:, idx[i]] + eps * sign[j] * Yu)
                else:
                    Xu = np.real(Orbit[0][:, idx[i]] - eps * sign[j] * Yu)

                stop_fun.direction = -(L - params[0])

                # Integrate unstable manifold forwards in time
                [SF, etf_u, states, times] = PCR3BP_propagator(Xu, et0, deltat,
                                                               prnt_out_dt, stop_fun, params[0], params[1], ang, L)

                states_u.append(states)
                times_u.append(times)
                SF_u = np.append([SF_u], [SF])

            if s:
                Ys = Orbit[2][:, 1]  # stable eigenvector

                Xs = np.real(Orbit[0][:, idx[i]] - eps * sign[j] * Ys)

                stop_fun.direction = (L - params[0])

                # Integrate stable manifold backwards in time
                [SF, etf, states, times] = PCR3BP_propagator(Xs, et0, -deltat,
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

    fig1, ax1 = plt.subplots()
    fig12, ax12 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig22, ax22 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig32, ax32 = plt.subplots()

    if len(states_s):
        for i in states_s:
            if abs(i[-1, 1]) < abs(pos[0] - mu1):
                if len(i[0, :0:-1]) % 2:
                    temp = fft(i[0, :1:-1] - pos[0] + i[1, :1:-1] * 1.j)
                    l = len(i[0, :1:-1])
                    signal_s.append(temp)
                else:
                    temp = fft(i[0, :0:-1] - pos[0] + i[1, :0:-1] * 1.j)
                    l = len(i[0, :0:-1])
                    signal_s.append(temp)
                xf = fftfreq(l, Data['prnt_out_dt'])
                xf = fftshift(xf)
                yf = fftshift(temp)
                ax1.plot(xf[abs(yf) > 0.01 * max(abs(yf))],
                         100.0 / l * np.real(yf[abs(yf) > 0.01 * max(abs(yf))]))
                ax12.plot(xf[abs(yf) > 0.01 * max(abs(yf))],
                          100.0 / l * np.imag(yf[abs(yf) > 0.01 * max(abs(yf))]))

    if len(states_u):
        for i in states_u:
            if abs(i[-1, 1]) < abs(pos[0] - mu1):
                if len(i[0, :-1]) % 2:
                    temp = fft(i[0, :-2] - pos[0] + i[1, :-2] * 1.j)
                    l = len(i[0, :-2])
                    signal_u.append(temp)
                else:
                    temp = fft(i[0, :-1] - pos[0] + i[1, :-1] * 1.j)
                    l = len(i[0, :-1])
                    signal_u.append(temp)
                xf = fftfreq(l, Data['prnt_out_dt'])
                xf = fftshift(xf)
                yf = fftshift(temp)
                ax2.plot(xf[abs(yf) > 0.01 * max(abs(yf))],
                         100.0 / l * np.real(yf[abs(yf) > 0.01 * max(abs(yf))]))
                ax22.plot(xf[abs(yf) > 0.01 * max(abs(yf))],
                          100.0 / l * np.imag(yf[abs(yf) > 0.01 * max(abs(yf))]))

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

            temp = fft(vecx - pos[0] + vecy * 1.j)
            l = len(vecx)

            signal_s.append(temp)

            xf = fftfreq(l, Data['prnt_out_dt'])
            xf = fftshift(xf)
            yf = fftshift(temp)
            ax3.plot(xf[abs(yf) > 0.01 * max(abs(yf))],
                     100.0 / l * np.real(yf[abs(yf) > 0.01 * max(abs(yf))]))
            ax32.plot(xf[abs(yf) > 0.01 * max(abs(yf))],
                      100.0 / l * np.imag(yf[abs(yf) > 0.01 * max(abs(yf))]))

    plt.style.use(h.THEME)
    # plt.show()


def plotm(mu1, mu2, pos, states_po, states_s, SF_s, states_u, SF_u, ang, angmin):
    fig1, ax1 = plt.subplots()
    ax1.plot(mu1, 0, 'ko')
    ax1.plot(pos[0], pos[1], 'ko')

    if len(states_po) == 2:
        ax1.plot(states_po[0][0][0], states_po[0][0][1], 'k')
        ax1.plot(states_po[1][0][0], states_po[1][0][1], 'k')
        ax1.plot(pos[2], pos[3], 'ko')
        size = [len(states_po[0][0][:, 0]), len(states_po[1][0][:, 0])]
    else:
        ax1.plot(states_po[0][0], states_po[0][1], 'k')
        size = [len(states_po[0][:, 0])]

    ax1.plot([mu1, mu1 + 1.25 * abs(pos[0] - mu1) * np.cos(angmin)],
             [0, 1.25 * abs(pos[0] - mu1) * np.sin(angmin)], 'g--')
    ax1.plot([mu1, mu1 + 1.25 * abs(pos[0] - mu1) * np.cos(angmin)],
             [0, -1.25 * abs(pos[0] - mu1) * np.sin(angmin)], 'g--')
    ax1.plot([mu1, mu1 + abs(pos[0] - mu1) * np.cos(ang * np.pi / 180)],
             [0, abs(pos[0] - mu1) * np.sin(ang * np.pi / 180)], 'r--')

    if max(size) == 4:

        fig2, ax2 = plt.subplots()
        for i_u in range(len(states_u)):
            ax2.plot(states_u[i_u][0], states_u[i_u][1], 'r')
        for i_s in range(len(states_s)):
            ax2.plot(states_s[i_s][0], states_s[i_s][1], 'b')
        ax2.plot(mu1, 0, 'ro')
        ax2.plot(pos[0], pos[1], 'ko')
        if len(states_po) == 2:
            ax1.plot(states_po[0][0][0], states_po[0][0][1], 'k')
            ax1.plot(states_po[1][0][0], states_po[1][0][1], 'k')
            ax1.plot(pos[2], pos[3], 'ko')
        else:
            ax1.plot(states_po[0][0], states_po[0][1], 'k')
        ax2.plot([mu1, mu1 + abs(pos[0] - mu1) * np.cos(ang * np.pi / 180)],
                 [0, abs(pos[0] - mu1) * np.sin(ang * np.pi / 180)], 'y--')

        h.log('The maximum tangential velocities registered in this Poincaré section are: ')

        u = 0
        s = 0
        if len(SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), :]) > 0:
            tangVels = max(abs(np.cos(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 2] \
                               + np.sin(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 3]))
            h.log('Max Vt (stable manifold):   %1.5f' % tangVels)

            xs = np.sqrt((SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 0] - mu1) ** 2 \
                         + SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 1] ** 2)
            ys = - np.sin(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 2] \
                 + np.cos(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 3]

            centers = [np.mean(xs), np.mean(ys)]
            atans = np.arctan2(ys - centers[1], xs - centers[0])
            keys = np.argsort(atans)

            s = 1

        if len(SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), :]) > 0:
            tangVelu = max(abs(np.cos(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 2] \
                               + np.sin(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 3]))
            h.log('Max Vt (unstable manifold): %1.5f' % tangVelu)

            xu = np.sqrt((SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 0] - mu1) ** 2 \
                         + SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 1] ** 2)
            yu = - np.sin(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 2] \
                 + np.cos(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 3]

            centeru = [np.mean(xu), np.mean(yu)]
            atanu = np.arctan2(yu - centeru[1], xu - centeru[0])
            keyu = np.argsort(atanu)

            u = 1

        if u or s:
            fig3, ax3 = plt.subplots()
            for j in range(len(SF_s)):
                if abs(SF_s[j, 1]) < abs(pos[0] - mu1):
                    ax3.plot(np.sqrt((SF_s[j, 0] - mu1) ** 2 + SF_s[j, 1] ** 2),
                             - np.sin(ang * np.pi / 180) * SF_s[j, 2] + np.cos(ang * np.pi / 180) * SF_s[j, 3],
                             'bo')
            for k in range(len(SF_u)):
                if abs(SF_u[k, 1]) < abs(pos[0] - mu1):
                    ax3.plot(np.sqrt((SF_u[k, 0] - mu1) ** 2 + SF_u[k, 1] ** 2),
                             - np.sin(ang * np.pi / 180) * SF_u[k, 2] + np.cos(ang * np.pi / 180) * SF_u[k, 3],
                             'ro')

            ax3.set(xlabel='Distance from 2nd primary')
            ax3.set(ylabel='Velocity normal to the plane')
            ax3.set(title=r'Poincaré section at angle %2.1f$^{\circ}$' % ang)

        if u:
            ax3.fill(xu[keyu], yu[keyu], facecolor='pink')
        if s:
            ax3.fill(xs[keys], ys[keys], facecolor='purple')

        plt.style.use(h.THEME)
        plt.show()

    else:

        fig2 = plt.figure()
        ax2 = Axes3D(fig2, auto_add_to_figure=False)
        fig2.add_axes(ax2)
        if size[0] == 4:
            for i_u in range(len(states_u)):
                ax2.plot(states_u[i_u][0], states_u[i_u][1], [0], 'r')
            for i_s in range(len(states_s)):
                ax2.plot(states_s[i_s][0], states_s[i_s][1], states_s[i_s][2], 'b')
            ax2.plot([mu1], [0], [0], 'ro')
            ax2.plot([pos[0]], [pos[1]], [0], 'ko')
            ax2.plot([pos[2]], [pos[3]], [0], 'ko')
            ax2.plot(states_po[0][0][0], states_po[0][0][1], [0], 'k')
            ax2.plot(states_po[1][0][0], states_po[1][0][1], states_po[1][0][2], 'k')

        elif size[-1] == 4:
            for i_u in range(len(states_u)):
                ax2.plot(states_u[i_u][0], states_u[i_u][1], states_u[i_u][2], 'r')
            for i_s in range(len(states_s)):
                ax2.plot(states_s[i_s][0], states_s[i_s][1], [0], 'b')
            ax2.plot([mu1], [0], [0], 'ro')
            ax2.plot([pos[0]], [pos[1]], [0], 'ko')
            ax2.plot([pos[2]], [pos[3]], [0], 'ko')
            ax2.plot(states_po[0][0][0], states_po[0][0][1], states_po[0][0][2], 'k')
            ax2.plot(states_po[1][0][0], states_po[1][0][1], [0], 'k')
        else:
            for i_u in range(len(states_u)):
                ax2.plot(states_u[i_u][0], states_u[i_u][1], states_u[i_u][2], 'r')
            for i_s in range(len(states_s)):
                ax2.plot(states_s[i_s][0], states_s[i_s][1], states_s[i_s][2], 'b')
            ax2.plot([mu1], [0], [0], 'ro')
            ax2.plot([pos[0]], [pos[1]], [0], 'ko')
            if len(states_po) == 2:
                ax2.plot(states_po[0][0][0], states_po[0][0][1], states_po[0][0][2], 'k')
                ax2.plot(states_po[1][0][0], states_po[1][0][1], states_po[1][0][2], 'k')
                ax2.plot([pos[2]], [pos[3]], [0], 'ko')
            else:
                ax2.plot(states_po[0][0], states_po[0][1], states_po[0][2], 'k')
        ax2.plot([mu1, mu1 + abs(pos[0] - mu1) * np.cos(ang * np.pi / 180)],
                 [0, abs(pos[0] - mu1) * np.sin(ang * np.pi / 180)],
                 [abs(pos[0] - mu1) / 2, abs(pos[0] - mu1) / 2], 'y--')
        ax2.plot([mu1, mu1 + abs(pos[0] - mu1) * np.cos(ang * np.pi / 180)],
                 [0, abs(pos[0] - mu1) * np.sin(ang * np.pi / 180)],
                 [-abs(pos[0] - mu1) / 2, -abs(pos[0] - mu1) / 2], 'y--')
        ax2.plot([mu1, mu1], [0, 0],
                 [-abs(pos[0] - mu1) / 2, abs(pos[0] - mu1) / 2], 'y--')
        ax2.plot([mu1 + abs(pos[0] - mu1) * np.cos(ang * np.pi / 180),
                  mu1 + abs(pos[0] - mu1) * np.cos(ang * np.pi / 180)],
                 [abs(pos[0] - mu1) * np.sin(ang * np.pi / 180),
                  abs(pos[0] - mu1) * np.sin(ang * np.pi / 180)],
                 [-abs(pos[0] - mu1) / 2, abs(pos[0] - mu1) / 2], 'y--')
        ax2.set(xlabel='x', ylabel='y', zlabel='z')

        h.log('The maximum tangential velocities registered in this Poincaré section are: ')

        u = 0
        s = 0

        if len(SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), :]) > 0:
            if size[-1] == 4:
                tangVels1 = max(abs(np.cos(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 2] \
                                    + np.sin(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 3]))
                h.log('Max Vt1 (stable manifold):   %1.5f' % tangVels1)
                tangVels2 = 0.00000
                h.log('Max Vt2 (stable manifold):   %1.5f' % tangVels2)

            else:
                tangVels1 = max(abs(np.cos(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 3] \
                                    + np.sin(ang * np.pi / 180) * SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 4]))
                h.log('Max Vt1 (stable manifold):   %1.5f' % tangVels1)
                tangVels2 = max(abs(SF_s[abs(SF_s[:, 1]) < abs(pos[0] - mu1), 5]))
                h.log('Max Vt2 (stable manifold):   %1.5f' % tangVels2)

            s = 1

        if len(SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), :]) > 0:
            if size[0] == 4:
                tangVelu1 = max(abs(np.cos(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 2] \
                                    + np.sin(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 3]))
                h.log('Max Vt1 (unstable manifold): %1.5f' % tangVelu1)
                tangVelu2 = 0.00000
                h.log('Max Vt2 (unstable manifold): %1.5f' % tangVelu2)

            else:
                tangVelu1 = max(abs(np.cos(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 3] \
                                    + np.sin(ang * np.pi / 180) * SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 4]))
                h.log('Max Vt1 (unstable manifold): %1.5f' % tangVelu1)
                tangVelu2 = max(abs(SF_u[abs(SF_u[:, 1]) < abs(pos[0] - mu1), 5]))
                h.log('Max Vt2 (unstable manifold): %1.5f' % tangVelu2)

            u = 1

        if u or s:
            fig3 = plt.figure()
            ax3 = Axes3D(fig3, auto_add_to_figure=False)
            fig3.add_axes(ax3)
            if size[0] == 4:
                for j in range(len(SF_s)):
                    if abs(SF_s[j, 1]) < abs(pos[0] - mu1):
                        ax3.plot([np.sqrt((SF_s[j, 0] - mu1) ** 2 + SF_s[j, 1] ** 2)],
                                 [SF_s[j, 2]],
                                 [-np.sin(ang * np.pi / 180) * SF_s[j, 3] + np.cos(ang * np.pi / 180) * SF_s[j, 4]],
                                 'bo')
                for k in range(len(SF_u)):
                    if abs(SF_u[k, 1]) < abs(pos[0] - mu1):
                        ax3.plot([np.sqrt((SF_u[k, 0] - mu1) ** 2 + SF_u[k, 1] ** 2)], [0],
                                 [-np.sin(ang * np.pi / 180) * SF_u[k, 2] + np.cos(ang * np.pi / 180) * SF_u[k, 3]],
                                 'ro')
            elif size[-1] == 4:
                for j in range(len(SF_s)):
                    if abs(SF_s[j, 1]) < abs(pos[0] - mu1):
                        ax3.plot([np.sqrt((SF_s[j, 0] - mu1) ** 2 + SF_s[j, 1] ** 2)], [0],
                                 [-np.sin(ang * np.pi / 180) * SF_s[j, 2] + np.cos(ang * np.pi / 180) * SF_s[j, 3]],
                                 'bo')
                for k in range(len(SF_u)):
                    if abs(SF_u[k, 1]) < abs(pos[0] - mu1):
                        ax3.plot([np.sqrt((SF_u[k, 0] - mu1) ** 2 + SF_u[k, 1] ** 2)],
                                 [SF_u[k, 2]],
                                 [-np.sin(ang * np.pi / 180) * SF_u[k, 3] + np.cos(ang * np.pi / 180) * SF_u[k, 4]],
                                 'ro')
            else:
                for j in range(len(SF_s)):
                    if abs(SF_s[j, 1]) < abs(pos[0] - mu1):
                        ax3.plot([np.sqrt((SF_s[j, 0] - mu1) ** 2 + SF_s[j, 1] ** 2)],
                                 [SF_s[j, 2]],
                                 [-np.sin(ang * np.pi / 180) * SF_s[j, 3] + np.cos(ang * np.pi / 180) * SF_s[j, 4]],
                                 'bo')
                for k in range(len(SF_u)):
                    if abs(SF_u[k, 1]) < abs(pos[0] - mu1):
                        ax3.plot([np.sqrt((SF_u[k, 0] - mu1) ** 2 + SF_u[k, 1] ** 2)],
                                 [SF_u[k, 2]],
                                 [-np.sin(ang * np.pi / 180) * SF_u[k, 3] + np.cos(ang * np.pi / 180) * SF_u[k, 4]],
                                 'ro')

            ax3.set(xlabel='Distance from 2nd primary (xy)')
            ax3.set(ylabel='Distance from 2nd primary (z)')
            ax3.set(zlabel='Velocity normal to the plane')
            ax3.set(title=r'Poincaré section at angle %2.1f$^{\circ}$' % ang)
        plt.style.use(h.THEME)
        plt.show()
