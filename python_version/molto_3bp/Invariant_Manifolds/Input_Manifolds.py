###########################################################################
# MANIFOLDS TEST CASES
# The example consists on a pair of manifolds intersecting around the L2 point
# in the Sun - Earth + Moon system.
###########################################################################
#
# Import functions required
#
from Manifolds import Manifolds
import numpy as np

# Main Functionalities (expected string)
# Type of orbit of interest, 'LY' == lyapunov orbit (2D)
#                            'HL' == halo orbit (3D)
type = 'LY'

# System analyzed
# Accepted values: 'SE', 'EM'
f = 'EM'

# Case choice (Orbital amplitude)
Ax       = 1.2e-3 # Lyapunov Orbit characterization

Az       = 100e3 # Halo Orbit characterization
phi      = 0
m        = 3    # 1 -> Northern variant of the orbit
                # 3 -> Southern variant of the orbit
LP       = 1    # 1 -> Lagrange point L1
                # 2 -> Lagrange point L2

poincSec = np.linspace(150, -150, 3)
                # Angle (in degrees) between the section required and
                # the +X semiplane taken from the 2nd primary

# Numerical parameters
npoints  = 10   # Number of points in the orbit to propagate manifolds from

d        = 1    # The program propagates the perturbation in the direction chosen
                # both directions 0, interior realm 1, exterior realm -1
                # Unexpected behaviour on halo orbits

branch   = 0    # Propagation of stable branch, unstable one or both
                # both branches 0, unstable branch -1, stable branch 1

prnt_out_dt = 0.001   # print time period

# In order to introduce data from the exterior, the program expects a .txt input

# Combining all input data into a single data structure
Input = {'type': type, 'mode': f, 'Ax_tgt': Ax, 'Az': Az, 'm': m, 'phi': phi,
    'LP': LP, 'poincSec': poincSec, 'npoints': npoints, 'd': d, 'branch': branch,
    'prnt_out_dt': prnt_out_dt}

Manifolds(Input)
