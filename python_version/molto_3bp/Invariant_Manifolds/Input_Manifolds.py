###########################################################################
# MANIFOLDS TEST CASES
# The example consists on a pair of manifolds intersecting around the L2 point
# in the Sun - Earth + Moon system.
###########################################################################
#
# Import functions required
#
from Main_Manifolds_SE import Main_Manifolds_SE

# Functionalities required by user
# Accepted values: True (same as any number) or False (same as 0)
# Any arrangment of values is valid but the program will ask the user for
# a source to find the data that the previous sections would have produced
f1  = True        # Initial guess generator and differential corrector
f2  = False       # Manifolds constructor and plotting tool
flags = [f1, f2]

# Case choice
Ax_tgt  = 1.4e-4

LP      = 2 # 1 -> Lagrange point L1
            # 2 -> Lagrange point L2

# Numerical parameters
npoints = 4          # Number of sections of the orbit to propagate manifolds

prnt_out_dt = 0.01   # print time period

# In order to introduce data externally, the program expects a .txt input

# Enter the .txt file on the current folder to retrieve the data
# A first line should contain the orbital period as the example given
# The format consists on 4 rows containing x, y, vx, vy from a Halo orbit
# Comments will not be processed according to the '#' Python standard if at the
# start of the line
IC = 'sample.txt'

# Combining all input data into a single data structure
Input = {'flags': flags, 'Ax_tgt': Ax_tgt, 'LP': LP, 'npoints': npoints,
    'prnt_out_dt': prnt_out_dt, 'IC': IC}

Main_Manifolds_SE(Input)
