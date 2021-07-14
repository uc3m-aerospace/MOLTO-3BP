import json

from Common.Controllers.OrbitController import OrbitController
import Utility.helpers as h
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
from scipy import linalg
import os
import Utility.exceptions as exe
import Common.Config.configurator as conf


class HaloController(OrbitController):
    HALO_ORBIT_SAMPLE_PATH = '\Externals\halo_orbit_sample.txt'

    def __init__(self, data):
        super().__init__(attributes={"data": json.dumps(data)})
        if not self._data:
            self._data = data

    def generator(self):
        '''
        HALO Orbit Generation
        :param self._data: the self._data passed into Halo Generator container as List
        :return:
        '''

        h.log('--- Halo Generator: 3rd order Richardson expansion ---')
        if self._data.get('mode') == 'SE':
            h.log('\nCR3BP: SUN-(EARTH+MOON) SYSTEM')
        else:
            h.log('\nCR3BP: EARTH-MOON SYSTEM')

        t = np.linspace(0, self._data.get('tf'), int(self._data.get('tf') / self._data.get('prnt_out_dt')) + 1)
        # s <-- Propagation time

        Tconversion = self._data.get('TSE') / (2 * np.pi)

        # Basic Variables Calculation

        rh = (self._data.get('mu') / 3) ** (1 / 3)
        gamma1ch = rh * (1 - 1 / 3 * rh - 1 / 9 * rh ** 2 - 1 / 27 * rh ** 3)
        gamma2ch = rh * (1 + 1 / 3 * rh - 1 / 9 * rh ** 2)

        # gamma polynomial
        poly = [1, -(3 - self._data.get('mu')), 3 - 2 * self._data.get('mu'), -self._data.get('mu'),
                2 * self._data.get('mu'),
                -self._data.get('mu')]  # for gamma 1 --> L1
        poly2 = [1, (3 - self._data.get('mu')), 3 - 2 * self._data.get('mu'), -self._data.get('mu'),
                 -2 * self._data.get('mu'),
                 -self._data.get('mu')]  # for gamma 2 --> L2

        root1 = np.roots(poly)
        root2 = np.roots(poly2)
        gamma1 = np.real(root1[np.imag(root1) == 0])
        gamma2 = np.real(root2[np.imag(root2) == 0])

        # cn coefficients calculation
        cn1 = lambda n: 1 / gamma1 ** 3 * (
                +1 ** n * self._data.get('mu') + (-1) ** n * ((1 - self._data.get('mu')) * gamma1 ** (n + 1)) / (
                1 - gamma1) ** (
                        n + 1))  # L1
        cn2 = lambda n: 1 / gamma2 ** 3 * (
                (-1) ** n * self._data.get('mu') + (-1) ** n * ((1 - self._data.get('mu')) * gamma2 ** (n + 1)) / (
                1 + gamma2) ** (
                        n + 1))  # L2

        if self._data.get('LP') == 1:
            c2 = cn1(2)
            c3 = cn1(3)
            c4 = cn1(4)
            xL = (1 - self._data.get('mu')) - gamma1  # Position of Lpoint w.r.t m1-m2 center of mass
            gammaC = gamma1
        elif self._data.get('LP') == 2:
            c2 = cn2(2)
            c3 = cn2(3)
            c4 = cn2(4)
            xL = (1 - self._data.get('mu')) + gamma2
            gammaC = gamma2
        else:
            config = conf.Configurator()
            if config.get('run.mode', 'debug').lower() != 'debug':
                raise exe.HaloOrbitLPError()
            else:
                exe.HaloOrbitLPError()

        r1 = self._data.get('Ltheory') * ((1 - self._data.get('mu')) - gamma1)
        r2 = self._data.get('Ltheory') * ((1 - self._data.get('mu')) + gamma2)

        # Eigenvalues Calculation

        wp = np.sqrt((2 - c2 + np.sqrt(9 * c2 ** 2 - 8 * c2)) / 2)
        wv = np.sqrt(c2)
        lmbda = wp  # For halo orbits

        # Variables - Linear solution

        kappa = (wp ** 2 + 1 + 2 * c2) / (2 * wp)
        kappa_C = 2 * lmbda / (lmbda ** 2 + 1 - c2)

        # Third-Order Richardson Expansion

        d1 = 3 * lmbda ** 2 / kappa * (kappa * (6 * lmbda ** 2 - 1) - 2 * lmbda)
        d2 = 8 * lmbda ** 2 / kappa * (kappa * (11 * lmbda ** 2 - 1) - 2 * lmbda)

        a21 = (3 * c3 * (kappa ** 2 - 2)) / (4 * (1 + 2 * c2))
        a22 = 3 * c3 / (4 * (1 + 2 * c2))
        a23 = -3 * c3 * lmbda / (4 * kappa * d1) * (3 * kappa ** 3 * lmbda - 6 * kappa * (kappa - lmbda) + 4)
        a24 = -3 * c3 * lmbda / (4 * kappa * d1) * (2 + 3 * kappa * lmbda)

        b21 = -3 * c3 * lmbda / (2 * d1) * (3 * kappa * lmbda - 4)
        b22 = 3 * c3 * lmbda / d1

        d21 = -c3 / (2 * lmbda ** 2)

        a31 = -9 * lmbda / (4 * d2) * (4 * c3 * (kappa * a23 - b21) + kappa * c4 * (4 + kappa ** 2)) + \
              (9 * lmbda ** 2 + 1 - c2) / (2 * d2) * (3 * c3 * (2 * a23 - kappa * b21) + c4 * (2 + 3 * kappa ** 2))
        a32 = -1 / d2 * (9 * lmbda / 4 * (4 * c3 * (kappa * a24 - b22) + kappa * c4) + \
                         3 / 2 * (9 * lmbda ** 2 + 1 - c2) * (c3 * (kappa * b22 + d21 - 2 * a24) - c4))

        b31 = 3 / (8 * d2) * (8 * lmbda * (3 * c3 * (kappa * b21 - 2 * a23) - c4 * (2 + 3 * kappa ** 2)) + \
                              (9 * lmbda ** 2 + 1 + 2 * c2) * (
                                      4 * c3 * (kappa * a23 - b21) + kappa * c4 * (4 + kappa ** 2)))
        b32 = 1 / d2 * (9 * lmbda * (c3 * (kappa * b22 + d21 - 2 * a24) - c4) + \
                        3 / 8 * (9 * lmbda ** 2 + 1 + 2 * c2) * (4 * c3 * (kappa * a24 - b22) + kappa * c4))

        d31 = 3 / (64 * lmbda ** 2) * (4 * c3 * a24 + c4)
        d32 = 3 / (64 * lmbda ** 2) * (4 * c3 * (a23 - d21) + c4 * (4 + kappa ** 2))

        s1 = (2 * lmbda * (lmbda * (1 + kappa ** 2) - 2 * kappa)) ** -1 * \
             (3 / 2 * c3 * (2 * a21 * (kappa ** 2 - 2) - a23 * (kappa ** 2 + 2) - 2 * kappa * b21) - 3 / 8 * c4 * (
                     3 * kappa ** 4 - 8 * kappa ** 2 + 8))
        s2 = (2 * lmbda * (lmbda * (1 + kappa ** 2) - 2 * kappa)) ** -1 * \
             (3 / 2 * c3 * (2 * a22 * (kappa ** 2 - 2) + a24 * (
                     kappa ** 2 + 2) + 2 * kappa * b22 + 5 * d21) + 3 / 8 * c4 * (12 - kappa ** 2))

        l1 = -3 / 2 * c3 * (2 * a21 + a23 + 5 * d21) - 3 / 8 * c4 * (12 - kappa ** 2) + 2 * lmbda ** 2 * s1
        l2 = 3 / 2 * c3 * (a24 - 2 * a22) + 9 / 8 * c4 + 2 * lmbda ** 2 * s2
        incre = wp ** 2 - wv ** 2

        # Relationships

        # Amplitude Constraint
        # Todo maybe normalise this
        Az = self._data.get('Az') / self._data.get('L')
        Ax = np.sqrt(-(incre + l2 * Az ** 2) / l1)  # km

        # Phase Angle Relationship
        if self._data.get('m') in [1, 3]:
            psi = self._data.get('phi') + self._data.get('m') * np.pi / 2  # rad
        else:
            config = conf.Configurator()
            if config.get('run.mode', 'debug') != 'debug':
                raise exe.HaloOrbitmError()
            else:
                exe.HaloOrbitmError()

        # Added variables - Lindstedt-Poncairé Method

        nu1 = 0
        nu2 = s1 * (Ax) ** 2 + s2 * (Az) ** 2
        nu = 1 + nu1 + nu2
        tau = nu * t

        tau1 = wp * tau + self._data.get('phi')
        deltam = 2 - self._data.get('m')

        T = 2 * np.pi / (wp * nu) * Tconversion  # Period Calculation (s)
        Tad = 2 * np.pi / (wp * nu)  # Period Calculation (-)
        Tdays = T / 3600 / 24  # s to days

        # EoM - 3rd Order Richardson Expansion

        xad = a21 * Ax ** 2 + a22 * Az ** 2 - Ax * np.cos(tau1) + \
              (a23 * Ax ** 2 - a24 * Az ** 2) * np.cos(2 * tau1) + (a31 * Ax ** 3 - a32 * Ax * Az ** 2) * np.cos(
            3 * tau1)
        yad = kappa * Ax * np.sin(tau1) + \
              (b21 * Ax ** 2 - b22 * Az ** 2) * np.sin(2 * tau1) + (b31 * Ax ** 3 - b32 * Ax * Az ** 2) * np.sin(
            3 * tau1)
        zad = deltam * Az * np.cos(tau1) + \
              deltam * d21 * Ax * Az * (np.cos(2 * tau1) - 3) + deltam * (d32 * Az * Ax ** 2 - d31 * Az ** 3) * np.cos(
            3 * tau1)

        # Introduce rx = f(xad)   ry = f(yad)   rz = f(zad)

        x = xad * gammaC + xL
        y = yad * gammaC
        z = zad * gammaC  # REVISAR

        # IC calculation --> For Initial Guess

        x0 = x[0]
        y0 = y[0]
        z0 = z[0]
        dt = t[1] - t[0]
        vx0 = (x[1] - x0) / dt
        vy0 = (y[1] - y0) / dt
        vz0 = (z[1] - z0) / dt

        # # Plots

        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)
        ax.plot(x, y, z)
        ax.plot(xL, 0, 0, 'ro')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        # 1 orbit display

        tdim = t * Tconversion
        fig2, (ax1, ax2, ax3) = plt.subplots(1, 3, constrained_layout=True)
        ax1.plot(x[tdim < T], y[tdim < T])
        ax1.plot(0, 0, 'bx')
        ax1.plot(xL, 0, 'ro')
        ax1.set(xlabel='x', ylabel='y')
        ax2.plot(x[tdim < T], z[tdim < T])
        ax2.plot(0, 0, 'bx')
        ax2.plot(xL, 0, 'ro')
        ax2.set(xlabel='x', ylabel='z')
        ax3.plot(y[tdim < T], z[tdim < T])
        ax3.plot(0, 0, 'bx')
        ax3.plot(0, 0, 'ro')
        ax3.set(xlabel='y', ylabel='z')
        if self._data.get('m') == 1:
            fig2.suptitle('Class I  Orbit  (L' + str(self._data.get('LP')) + ')')
        else:
            fig2.suptitle('Class II  Orbit  (L' + str(self._data.get('LP')) + ')')
        plt.show()

        # # Displays
        h.log('-- Results ---')
        h.log('Ax = %.1f' % (Ax * self._data.get('L')) + ' km   Az = ' + str(self._data.get('Az')) + ' km')
        h.log('T (days) = %.5f' % Tdays)
        if self._data.get('m') == 1:
            h.log('Case: Northern (L' + str(self._data.get('LP')) + ')')
        else:
            h.log('Case: Southern (L' + str(self._data.get('LP')) + ')')

        # IC guess display
        h.log('\n--- Initial Guess Display ---')
        h.log('x0  = %.20f' % x0)
        h.log('y0  = %.20f' % y0)
        h.log('z0  = %.20f' % z0)
        h.log('vx0 = %.20f' % vx0)
        h.log('vy0 = %.20f' % vy0)
        h.log('vz0 = %.20f\n' % vz0)

        # Time guesses
        h.log('--- Time Guess ---')
        h.log('T  = %.20f' % Tad)
        h.log('T/2  = %.20f\n' % (Tad / 2))

        if self._data.get('flags')[1]:
            self._data.get['IC'] = np.array([x0, z0, vy0])
            (states_po, T_po, eigvec) = self.num_comp(self._data)

            return (states_po, T_po, eigvec)
        else:

            text = '# self._data Produced by Halo Generator #\n' + '# mode = ' + str(self._data.get('mode')) + \
                   ' LP = ' + str(self._data.get('LP')) + ' m = ' + str(self._data.get('m')) + ' phi = ' + \
                   str(self._data.get('phi')) + ' Az = ' + str(self._data.get('Az')) + '\n' + \
                   'x0  = %.20f\n' % x0 + 'y0  = %.20f\n' % y0 + 'z0  = %.20f\n' % z0 + \
                   'vx0 = %.20f\n' % vx0 + 'vy0 = %.20f\n' % vy0 + 'vz0 = %.20f\n' % vz0
            fid = open(self.HALO_ORBIT_SAMPLE_PATH, 'w')
            fid.write(text)
            fid.close()

    def plot(self, plot_1=True, plot_2=True) -> (Any, Any, Any):
        ##### HALO ORBITS PLOTTING TOOL #####
        #
        # Importing required functions
        #

        # Initial conditions
        if self._data.get('flags')[0] or self._data.get('flags')[1] or self._data.get('method') == 'insitu':
            [x0, z0, vy0] = self._data.get('IC')
        else:
            fid = open(self._data.get('IC'), 'r')
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
                config = conf.Configurator()
                if config.get('run.mode', 'debug').lower() != 'debug':
                    raise exe.HaloNumCompError()
                else:
                    exe.HaloNumCompError()

        q0 = np.array([x0, 0, z0, 0, vy0, 0])
        tspan = np.linspace(0, self._data.get('tf'), int(self._data.get('tf') / self._data.get('prnt_out_dt')) + 1)

        sol = solve_ivp(h.three_body_problem, [0, self._data.get('tf')],
                        q0, t_eval=tspan, args=(self._data.get('mu'),),
                        atol=1e-15,
                        rtol=1e-10)
        times_po = sol.t
        states_po = sol.y
        x = states_po[0]
        y = states_po[1]
        z = states_po[2]

        q0 = np.zeros(42)
        q0[:6] = [x0, 0, z0, 0, vy0, 0]
        phi0 = np.identity(6)
        q0[6:] = phi0.ravel()

        sol = solve_ivp(h.diff_correction, [0, self._data.get('tf')],
                        q0, args=(self._data.get('mu'),),
                        atol=1e-15,
                        rtol=1e-10)
        q = sol.y
        Phi = q[6:, -1].reshape(6, -1)

        [eigval, eigvec] = linalg.eig(Phi)

        # Figures

        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=True)
        fig.add_axes(ax)
        ax.plot(x, y, z)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        fig2, (ax1, ax2, ax3) = plt.subplots(1, 3, constrained_layout=True)
        ax1.plot(x, y)
        ax1.set(xlabel='x', ylabel='y')
        ax2.plot(x, z)
        ax2.set(xlabel='x', ylabel='z')
        ax3.plot(y, z)
        ax3.set(xlabel='y', ylabel='z')
        plt.show()

        return (states_po, self._data.get('tf'), eigvec[:, np.imag(eigval) == 0])

    def num_comp(self):
        '''
        #### HALO ORBITS NUMERICAL COMPUTATION #####

        Importing required functions


        IC guess

        Use Results from Halo_Generator.py or sample.txt

        --- Initial Guess ---
        '''

        if self._data.get('flags')[0] or self._data.get('method') == 'insitu':
            [x0, z0, vy0] = self._data.get('IC')
        else:
            fid = open(self._data.get('IC'), 'r')
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
                config = conf.Configurator()
                if config.get('run.mode', 'debug').lower() != 'debug':
                    raise exe.HaloNumCompError()
                else:
                    exe.HaloNumCompError()

        # Recall y0, vx0, vz0 = 0

        ## Differential Correction
        for i in range(1, self._data.get('nmax') + 1):

            # IC vector preparation

            # x' = f(x)    IC
            q0 = np.zeros(42)  # Initial Conditions
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
            sol = solve_ivp(h.diff_correction, [0, self._data.get('tf')],
                            q0, events=EventHalo, args=(self._data.get('mu'),),
                            atol=1e-8,
                            rtol=1e-5)
            q = sol.y
            t = sol.t

            # Extracting solution
            xfvec = q[:6, -1]
            xbfvec = q[:6, -2]

            Phifvec = q[6:, -1]
            Phifmat = Phifvec.reshape(6, -1)  # to matrix form

            # Desired values at tf: vxf,vzf = 0   (yf=0 already with stop event)
            # Desired values             % Final values
            vxfdes = 0
            vzfdes = 0
            vxf = xfvec[3]
            vzf = xfvec[5]

            ydot = xfvec[4]
            xdotdot = (vxf - xbfvec[3]) / (t[-1] - t[-2])
            zdotdot = (vzf - xbfvec[5]) / (t[-1] - t[-2])

            # Delta x
            dvx = vxfdes - vxf
            dvz = vzfdes - vzf
            B = np.array([[dvx], [dvz]])
            D = np.array([[xdotdot], [zdotdot]])
            E = np.array([Phifmat[1, 0], Phifmat[1, 4]])

            # Check of IC

            err1 = abs(dvx)
            err2 = abs(dvz)

            if (err1 <= self._data.get('tol')) and (err2 <= self._data.get('tol')):
                break
            else:

                # Update IC --- Ax=B
                A = np.array(Phifmat[np.array([3, 3, 5, 5]),
                                     np.array([0, 4, 0, 4])].reshape((2, 2)))
                C = A - 1 / ydot * D * E
                dxvec0 = linalg.solve(C, B)  # Solve inverting C
                dx0 = dxvec0[0]
                dvy0 = dxvec0[1]

                x0 = x0 + dx0
                vy0 = vy0 + dvy0

        ## Solution

        h.log('--- Halo Generator: Numerical Computation ---\n')
        if (err1 <= self._data.get('tol')) and (err2 <= self._data.get('tol')):
            h.log('Nº of iterations to converge: ' + str(i))
            h.log('\n--- Solution ---')
            h.log('x0  = %.20f' % x0)
            h.log('y0  = 0.00')
            h.log('z0  = %.20f' % z0)
            h.log('vx0 = 0.00')
            h.log('vy0 = %.20f' % vy0)
            h.log('vz0 = 0.00\n')
            h.log('--- Orbit Period ---')
            h.log('T = %.20f' % (t[-1] * 2))
            h.log('T/2 = %.20f\n' % t[-1])
        else:
            h.log('The program has not converged!')
            h.log('err1  = ' + str(err1))
            h.log('err2  = ' + str(err2))
            h.log('\nNº of iterations done: ' + str(i))
            h.log('Tolerance: ' + str(self._data.get('tol')))
            h.log('Try modifying the initial guess IC ...')
            h.log('  ...or modifying the number of iterations and/or the tolerance')

        if self._data.get('flags')[2]:
            self._data['IC'] = np.array([x0[0], z0, vy0[0]])
            self._data['tf'] = t[-1] * 2

            (states_po, T_po, eigvec) = self.plot()

            return (states_po, T_po, eigvec)
        else:
            text = '# Data Produced by Halo Numerical Computation #\n' + '# opt = ' + \
                   str(self._data.get('opt')) + ' LP = ' + str(self._data.get('LP')) + ' m = ' + \
                   str(self._data.get('m')) + ' phi = ' + str(self._data.get('phi')) + ' Az = ' + \
                   str(self._data.get('Az')) + '\n' + 'x0  = %.20f\n' % x0 + 'z0  = %.20f\n' % z0 + \
                   'vy0 = %.20f\n' % vy0
            fid = open(self.HALO_ORBIT_SAMPLE_PATH, 'w')
            fid.write(text)
            fid.close()

    def main(self):

        if (type(self._data['flags'][0]) is not bool and type(self._data['flags'][0]) is not int) \
                or (type(self._data['flags'][1]) is not bool and type(self._data['flags'][1]) is not int) \
                or (type(self._data['flags'][2]) is not bool and type(self._data['flags'][2]) is not int) \
                or (
                self._data['flags'][0] == self._data['flags'][1] and self._data['flags'][0] == self._data['flags'][2] \
                and self._data['flags'][0] == 0):
            config = conf.Configurator()
            if config.get('run.mode', 'debug').lower() != 'debug':
                raise exe.HaloMainFlagsError()
            else:
                exe.HaloMainFlagsError()

        if self._data['flags'][0]:
            if sum(self._data['flags']) == 3:
                (states_po, T_po, eigvec) = \
                    self.generator()
                return (states_po, T_po, eigvec)
            else:
                self.generator()
        else:
            if self._data['method'] != 'insitu' and self._data['method'] != 'text':
                config = conf.Configurator()
                if config.get('run.mode', 'debug').lower() != 'debug':
                    raise exe.HaloMainMethodError()
                else:
                    exe.HaloMainMethodError()

            if self._data['flags'][1]:
                self.num_comp()
            elif self._data['flags'][2]:
                self.plot()
