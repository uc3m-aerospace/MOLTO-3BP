from spiceypy import spice_error_check
import Utility.decorators as d
import Utility.exceptions as ex
import Common.Config.configurator as c
import numpy as np
import spiceypy as spice
import os
from Common.Config.configurator import Configurator as c
CONFIG = 'debug'

THEME = c().get('plots.theme')

@d.retry_and_log(retry_attempts=1)
def load_kernel(kernel="gm_de431.tpc", raised=True):
    '''
    Load specific Kernel specified as default to gm_de431
    Beware the kernels should be loaded in the Externals/Kernels folder
    :return: True if loaded else will throw error
    '''


    workdir = os.getcwd()

    KERNEL_DIR = workdir + '/Externals/Kernels/'

    try:
        spice.furnsh(KERNEL_DIR + kernel)
    except Exception as e:
        log(repr(e))
        config = c.Configurator()
        if raised and config.get('run.mode', 'DEBUG') != 'DEBUG':
            raise ex.UnknownKernelType()
        else:
            ex.UnknownKernelType()

        return False

    return True


def log(value, is_for_slack=False, is_critical=False):
    print(value)
    value = str(value).replace("'", "")


def three_body_problem(t, q, mu):
    '''
    For the Halo computation solving the Diff Eq 3BP
    :param t:
    :param q:
    :param mu:
    :return:
    '''
    #### State Transition Matrix for CR3BP Computation ####
    # This function will solve the following ODEs system
    #                x'  = f(x)                   6x1

    import numpy as np

    # State vector

    # Dyn var

    x = q[0]
    y = q[1]
    z = q[2]
    vx = q[3]
    vy = q[4]
    vz = q[5]

    # Data

    mu1 = 1 - mu
    mu2 = mu

    ## Matrix Computation

    # # x' = f(x)

    Ux = mu2 * (mu1 - x) - mu1 * (mu2 + x) \
         - (mu2 * (mu1 - x)) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
         + (mu1 * (mu2 + x)) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (3 / 2)
    Uy = (mu1 * y) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) - mu2 * y - mu1 * y \
         + (mu2 * y) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (3 / 2)
    Uz = (mu1 * z) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
         + (mu2 * z) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (3 / 2)

    ## Ode Work

    # x' = f(x)

    xdot = vx
    ydot = vy
    zdot = vz
    xdotdot = 2 * vy - Ux
    ydotdot = -2 * vx - Uy
    zdotdot = -Uz

    return [xdot, ydot, zdot, xdotdot, ydotdot, zdotdot]


def diff_correction(t, q, mu):
    '''
    Differential correction/optimisation for Halo problem
    :param t:
    :param q:
    :param mu:
    :return:
    '''
    #### State Transition Matrix for CR3BP Computation ####

    # This function will solve the following ODEs system
    #                x'  = f(x)                   6x1
    #               Phi' = Df(x)*Phi(t,t0)        6x6

    x = q[0]
    y = q[1]
    z = q[2]

    # Derivative vector preallocation
    dqdt = np.zeros(len(q))

    # Dyn var
    dqdt[:6] = three_body_problem(t, q[:6], mu)

    # STM (State Transition Matrix)
    Vec = q[6:]  # Aux vector
    Phi = Vec.reshape((6, -1))  # Matrix initialized and assigned

    mu1 = 1 - mu
    mu2 = mu

    ## Matrix Computation

    # # Phi' = Df(x)*Phi(t,t0)

    # Matrix U: Uxx, Uyy, Uzz .... (Sign - is still involved in this matrix)
    # U is equivalent to -U matrix from Koon ref. (P161)

    m11 = (3 * mu1 * (2 * mu2 + 2 * x) ** 2) / (4 * ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)) \
          - mu1 / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
          - mu2 / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
          + (3 * mu2 * (2 * mu1 - 2 * x) ** 2) / (4 * ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)) + 1

    m12 = (3 * mu1 * y * (2 * mu2 + 2 * x)) / (2 * ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)) \
          - (3 * mu2 * y * (2 * mu1 - 2 * x)) / (2 * ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2))

    m13 = (3 * mu1 * z * (2 * mu2 + 2 * x)) / (2 * ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)) \
          - (3 * mu2 * z * (2 * mu1 - 2 * x)) / (2 * ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2))

    m21 = (3 * mu1 * y * (2 * mu2 + 2 * x)) / (2 * ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)) \
          - (3 * mu2 * y * (2 * mu1 - 2 * x)) / (2 * ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2))

    m22 = (3 * mu2 * y ** 2) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2) \
          - mu1 / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
          - mu2 / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
          + (3 * mu1 * y ** 2) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2) + 1

    m23 = (3 * mu2 * y * z) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2) \
          + (3 * mu1 * y * z) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)

    m31 = (3 * mu1 * z * (2 * mu2 + 2 * x)) / (2 * ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)) \
          - (3 * mu2 * z * (2 * mu1 - 2 * x)) / (2 * ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2))

    m32 = (3 * mu2 * y * z) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2) \
          + (3 * mu1 * y * z) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)

    m33 = (3 * mu2 * z ** 2) / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (5 / 2) \
          - mu1 / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
          - mu2 / ((mu1 - x) ** 2 + y ** 2 + z ** 2) ** (3 / 2) \
          + (3 * mu1 * z ** 2) / ((mu2 + x) ** 2 + y ** 2 + z ** 2) ** (5 / 2)

    U = np.array([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

    # Omega Matrix + Identity matrix (I3) + zeros 3x3 matrix (Z3)
    # Equivalent to 2*Omega matrix in Koon ref.(P161)

    Omega = np.array([[0, 2, 0], [-2, 0, 0], [0, 0, 0]])

    I3 = np.identity(3)
    Z3 = np.zeros((3, 3))

    # Df(x) Matrix
    Df = np.append(np.append(Z3, I3, axis=1),  # 6x6 Matrix
                   np.append(U, Omega, axis=1), axis=0)

    ## Ode Work

    # Phi' = Df(x)*Phi(t,t0)

    Phidot = Df @ Phi

    ## State Vector Derivative

    PhidotAux = Phidot.ravel()

    dqdt[6:] = PhidotAux

    return dqdt


def poinc_crossing(et, state, mu1, mu2, ang, L):
    '''
    Crossing Determinant
    :param et:
    :param state:
    :param mu1:
    :param mu2:
    :param ang:
    :param L:
    :return:
    '''

    import numpy as np

    [x, y] = state[:2]

    d = ((L - mu1)*np.sin(ang*np.pi/180)*(x-mu1)\
        - (L - mu1)*np.cos(ang*np.pi/180)*y)

    xprime = np.cos(ang*np.pi/180)*(x-mu1) + np.sin(ang*np.pi/180)*y

    cond = xprime > 0

    if abs(y) > abs(L-mu1):
        return 0
    elif cond:
        return d
    else:
        if ang - 180 > 0:
            return -mu1
        else:
            return mu1
