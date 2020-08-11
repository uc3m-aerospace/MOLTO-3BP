###########################################################################
# MANIFOLDS TEST CASES
# The example consists on a pair of manifolds intersecting around the L2 point
# in the Sun - Earth + Moon system.
###########################################################################
#
# Import functions required
#
from Manifolds import Manifolds

# Main Functionalities (expected string)
# Type of orbit of interest, 'LY' == lyapunov orbit (2D)
#                            'HL' == halo orbit (3D)
type = 'HL'

# System analyzed
# Accepted values: 'SE', 'EM'
f = 'EM'

# Case choice (Orbital amplitude)
Ax       = 1.2e-3 # Lyapunov Orbit characterization

Az       = 95e3 # Halo
m        = 3    # Orbit
phi      = 0    # characterization

LP       = 1    # 1 -> Lagrange point L1
                # 2 -> Lagrange point L2

poincSec = 0     # Angle (in degrees) between the section required and
                 # the +X semiplane taken from the 2nd primary

# Numerical parameters
npoints  = 10         # Number of points in the orbit to propagate manifolds from
                      # The program propagates the perturbation in both directions

prnt_out_dt = 0.001   # print time period

# In order to introduce data from the exterior, the program expects a .txt input

# Combining all input data into a single data structure
Input = {'type': type, 'mode': f, 'Ax_tgt': Ax, 'Az': Az, 'm': m, 'phi': phi,
    'LP': LP, 'poincSec': poincSec, 'npoints': npoints,
    'prnt_out_dt': prnt_out_dt}

Manifolds(Input)
