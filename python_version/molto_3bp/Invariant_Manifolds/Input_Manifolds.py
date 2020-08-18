###########################################################################
# MANIFOLDS TEST CASES
# Use this sheet as format for the study of a single type of orbit
# For sequential analysis proceed to set input_seq = True
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
Ax       = 1.15e-3 # Lyapunov Orbit characterization

Az       = 200e3 # Halo Orbit characterization
phi      = 0
m        = 3    # 1 -> Northern variant of the orbit
                # 3 -> Southern variant of the orbit
LP       = 1    # 1 -> Lagrange point L1
                # 2 -> Lagrange point L2

poincSec = np.array([0])
                # Angle (in degrees) between the section required and
                # the +X semiplane taken from the 2nd primary

# Numerical parameters
npoints  = 10   # Number of points in the orbit to propagate manifolds from

d        = 0    # The program propagates the perturbation in the direction chosen
                # both directions 0, interior realm 1, exterior realm -1
                # Unexpected behaviour on halo orbits

branch   = 0    # Propagation of stable branch, unstable one or both
                # both branches 0, unstable branch -1, stable branch 1

prnt_out_dt = 0.001   # print time period

input_seq = 0   # Sets the program to evaluate a mission from orbit 1 .. n
                # Following intersections in 1 .. n Poincaré sections
                # The only values overwritten by this parameter are the case
                # choice variables and d and branch parameters to create the chain

# In order to introduce data from the exterior, the program expects a .txt input

# Combining all input data into a single data structure
Input = {'type': type, 'mode': f, 'Ax_tgt': Ax, 'Az': Az, 'm': m, 'phi': phi,
    'LP': LP, 'poincSec': poincSec, 'npoints': npoints, 'd': d, 'branch': branch,
    'prnt_out_dt': prnt_out_dt, 'input_seq': input_seq}

Manifolds(Input)
