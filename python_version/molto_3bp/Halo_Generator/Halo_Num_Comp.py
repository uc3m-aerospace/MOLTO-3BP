##### HALO ORBITS NUMERICAL COMPUTATION #####
#
# Importing required functions
#
import numpy as np
from scipy.integrate import solve_ivp
from scipy import linalg
from intFun import DiffCorrection

# Numerical Computation characteristics
nmax = 50 # Max. number of iteration
tol  = 1e-15 # Tolerance error

# Orbit guess (Establish max for ode work)
tf = 3

## IC guess

# Use Results from Halo_Generator.py -> IC.txt

#--- Initial Guess Display ---
x0  = 1.1124550077766104
y0  = 0
z0  = 0.035680331960522345
vx0 = 0.0001677345614018
vy0 = 0.20156708661850475
vz0 = -0.0010217302462787591

# Recall!!! y0,vx0,vz0 = 0 !!!!!
y0 = 0; vx0 = 0; vz0 = 0

## Differential Correction
for i in range(1, nmax+1):

# IC vector preparation

# x' = f(x)    IC
    q0 = np.zeros(42) # Initial Conditions
    q0[:6] = [x0, y0, z0, vx0, vy0, vz0]

# Phi' = Df(x)*Phi(t,t0)     IC
    Phi0 = np.identity(6)
    q0[6:] = Phi0.ravel()

# Ode45 - State Transition Matrix Computation

# Event to stop ode45, x-z plane cross (at T/2)
    def EventHalo(t, q):
        if t > 0.:
            return q[1]
        else:
            return 1
    EventHalo.terminal = True
    sol = solve_ivp(DiffCorrection, [0, tf],
        q0, events = EventHalo,
        atol = 1e-9,
        rtol = 1e-6)
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

    if (err1<=tol) and (err2<=tol):
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
if (err1<=tol) and (err2<=tol):
    print('Nº of iterations to converge: ' + str(i))
    print('\n--- Solution ---')
    print('x0  = %.20f;' % x0)
    print('y0  = %.20f;' % y0)
    print('z0  = %.20f;' % z0)
    print('vx0  = %.20f;' % vx0)
    print('vy0  = %.20f;' % vy0)
    print('vz0  = %.20f;\n' % vz0)
    print('--- Orbit Period ---')
    print('T/2 = %.20f;' % t[-1])
    print('T = %.20f;' % (t[-1]*2))
else:
    print('The program has not converged!')
    print('err1  = ' + str(err1))
    print('err2  = ' + str(err2))
    print('\nNº of iterations done: ' + str(i))
    print('Tolerance: ' + str(tol))
    print('Try modifying the initial guess IC ...')
    print('  ...or modifying the number of iterations and/or the tolerance')
