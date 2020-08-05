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
# Accepted values: 'SE', 'EM', 'ALL'
# In 'ALL' mode the program will set both flags to True since in order to
# compose the mission it is necessary to build both manifolds and join them
f = 'SE'

# Specific Functionalities required by user
# Accepted values: True (same as any number) or False (same as 0)
# Any arrangment of values is valid but the program will ask the user for
# a source to find the data that the previous sections would have produced
f1  = 0        # Initial guess generator, differential corrector and orbit propagator
f2  = 1        # Manifolds constructor and plotting tool
flags = [f1, f2]

# Case choice
Ax_tgt   = 2e-4

LP       = 1    # 1 -> Lagrange point L1
                # 2 -> Lagrange point L2

poincSec = 0     # Angle (in degrees) between the section required and
                 # the +X semiplane taken from the 2nd primary

# Numerical parameters
npoints  = 10         # Number of points in the orbit to propagate manifolds from
                      # The program propagates the perturbation in both directions

prnt_out_dt = 0.01    # print time period (if f1 = False, must be consistent
                      # with the sampling rate of the input data!)

# In order to introduce data from the exterior, the program expects a .txt input

# Enter the .txt file on the current folder to retrieve the data
# A first line should contain the orbital period as the example given
# The format consists on 4 rows containing x, y, vx, vy from a Halo orbit
# Comments will not be processed according to the '#' Python standard if at the
# start of the line
IC = 'sample.txt'

# Combining all input data into a single data structure
Input = {'mode': f, 'flags': flags, 'Ax_tgt': Ax_tgt, 'LP': LP, 'poincSec': poincSec,
    'npoints': npoints, 'prnt_out_dt': prnt_out_dt, 'IC': IC}

Manifolds(Input)
