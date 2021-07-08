def ThreeBodyProp(t, q, mu):
#### State Transition Matrix for CR3BP Computation ####
# This function will solve the following ODEs system
#                x'  = f(x)                   6x1

    import numpy as np

# State vector

# Dyn var

    x  = q[0]
    y  = q[1]
    z  = q[2]
    vx = q[3]
    vy = q[4]
    vz = q[5]

# Data

    mu1 = 1 - mu
    mu2 = mu

## Matrix Computation

# # x' = f(x)

    Ux = mu2*(mu1 - x) - mu1*(mu2 + x)\
        - (mu2*(mu1 - x))/((mu1 - x)**2 + y**2 + z**2)**(3/2)\
        + (mu1*(mu2 + x))/((mu2 + x)**2 + y**2 + z**2)**(3/2)
    Uy = (mu1*y)/((mu2 + x)**2 + y**2 + z**2)**(3/2) - mu2*y - mu1*y\
        + (mu2*y)/((mu1 - x)**2 + y**2 + z**2)**(3/2)
    Uz = (mu1*z)/((mu2 + x)**2 + y**2 + z**2)**(3/2)\
        + (mu2*z)/((mu1 - x)**2 + y**2 + z**2)**(3/2)

## Ode Work

# x' = f(x)

    xdot = vx
    ydot = vy
    zdot = vz
    xdotdot =  2*vy-Ux
    ydotdot = -2*vx-Uy
    zdotdot = -Uz

    return [xdot, ydot, zdot, xdotdot, ydotdot, zdotdot]

def DiffCorrection(t, q, mu):
#### State Transition Matrix for CR3BP Computation ####

# This function will solve the following ODEs system
#                x'  = f(x)                   6x1
#               Phi' = Df(x)*Phi(t,t0)        6x6
    import numpy as np

    x  = q[0]
    y  = q[1]
    z  = q[2]

# Derivative vector preallocation
    dqdt = np.zeros(len(q))

# Dyn var
    dqdt[:6] = ThreeBodyProp(t, q[:6], mu)

# STM (State Transition Matrix)
    Vec = q[6:] # Aux vector
    Phi = Vec.reshape((6,-1)) # Matrix initialized and assigned

    mu1 = 1 - mu
    mu2 = mu

## Matrix Computation

# # Phi' = Df(x)*Phi(t,t0)

# Matrix U: Uxx, Uyy, Uzz .... (Sign - is still involved in this matrix)
# U is equivalent to -U matrix from Koon ref. (P161)

    m11 = (3*mu1*(2*mu2 + 2*x)**2)/(4*((mu2 + x)**2 + y**2 + z**2)**(5/2))\
        - mu1/((mu2 + x)**2 + y**2 + z**2)**(3/2)\
        - mu2/((mu1 - x)**2 + y**2 + z**2)**(3/2)\
        + (3*mu2*(2*mu1 - 2*x)**2)/(4*((mu1 - x)**2 + y**2 + z**2)**(5/2)) + 1

    m12 = (3*mu1*y*(2*mu2 + 2*x))/(2*((mu2 + x)**2 + y**2 + z**2)**(5/2))\
        - (3*mu2*y*(2*mu1 - 2*x))/(2*((mu1 - x)**2 + y**2 + z**2)**(5/2))

    m13 = (3*mu1*z*(2*mu2 + 2*x))/(2*((mu2 + x)**2 + y**2 + z**2)**(5/2))\
        - (3*mu2*z*(2*mu1 - 2*x))/(2*((mu1 - x)**2 + y**2 + z**2)**(5/2))

    m21 =(3*mu1*y*(2*mu2 + 2*x))/(2*((mu2 + x)**2 + y**2 + z**2)**(5/2))\
        - (3*mu2*y*(2*mu1 - 2*x))/(2*((mu1 - x)**2 + y**2 + z**2)**(5/2))

    m22 = (3*mu2*y**2)/((mu1 - x)**2 + y**2 + z**2)**(5/2)\
        - mu1/((mu2 + x)**2 + y**2 + z**2)**(3/2)\
        - mu2/((mu1 - x)**2 + y**2 + z**2)**(3/2)\
        + (3*mu1*y**2)/((mu2 + x)**2 + y**2 + z**2)**(5/2) + 1

    m23 = (3*mu2*y*z)/((mu1 - x)**2 + y**2 + z**2)**(5/2)\
        + (3*mu1*y*z)/((mu2 + x)**2 + y**2 + z**2)**(5/2)

    m31 = (3*mu1*z*(2*mu2 + 2*x))/(2*((mu2 + x)**2 + y**2 + z**2)**(5/2))\
        - (3*mu2*z*(2*mu1 - 2*x))/(2*((mu1 - x)**2 + y**2 + z**2)**(5/2))

    m32 = (3*mu2*y*z)/((mu1 - x)**2 + y**2 + z**2)**(5/2)\
        + (3*mu1*y*z)/((mu2 + x)**2 + y**2 + z**2)**(5/2)

    m33 = (3*mu2*z**2)/((mu1 - x)**2 + y**2 + z**2)**(5/2)\
        - mu1/((mu2 + x)**2 + y**2 + z**2)**(3/2)\
        - mu2/((mu1 - x)**2 + y**2 + z**2)**(3/2)\
        + (3*mu1*z**2)/((mu2 + x)**2 + y**2 + z**2)**(5/2)

    U = np.array([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

# Omega Matrix + Identity matrix (I3) + zeros 3x3 matrix (Z3)
# Equivalent to 2*Omega matrix in Koon ref.(P161)

    Omega = np.array([[0, 2, 0], [-2, 0, 0], [0, 0, 0]])

    I3 = np.identity(3); Z3 = np.zeros((3,3))

# Df(x) Matrix
    Df = np.append(np.append(Z3, I3, axis = 1), # 6x6 Matrix
        np.append(U, Omega, axis = 1), axis = 0)

## Ode Work

# Phi' = Df(x)*Phi(t,t0)

    Phidot = Df @ Phi

## State Vector Derivative

    PhidotAux = Phidot.ravel()

    dqdt[6:] = PhidotAux

    return dqdt
