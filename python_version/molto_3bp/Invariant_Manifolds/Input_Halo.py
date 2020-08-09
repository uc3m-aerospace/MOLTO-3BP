###########################################################################
# HALO TEST CASES
# The example consists on a Halo orbit around the L1 point in the Sun -
# Earth + Moon system over the north. (phi = 0; Az = 95000 km)
###########################################################################
#
# Import functions required
#
import numpy as np
from Halo_Orbits.Halo_Main import Halo_Main

# Functionalities required by user
# Accepted values: True (same as any number) or False (same as 0)
# Any arrangment of values is valid but the program will ask the user for
# a source to find the data that the previous sections would have produced
f1  = 1      # Initial guess generator from 3rd Order Richardson Expansion
f2  = 1      # Differential Correction module to refine Initial guess
f3  = 1      # Propagator and plotting tool
flags = [f1, f2, f3]

# Time period of study
tf     = 4      # Fix a value close to the expected period of the orbit
                # in order to improve performance
nsteps = 1e5

# Case choice
opt = 1 # This first one is required for any use of the program
        # Valid choices: 1 -> Sun - Earth+Moon System
        #                2 -> Earth - Moon Syst
LP  = 1 # 1 -> Lagrange point L1
        # 2 -> Lagrange point L2
m   = 1 # 1 -> Northern variant of the orbit
        # 3 -> Southern variant of the orbit

# Initial Orbital parameters
phi = 0         # rad
Az  = 95e3      # km

# Numerical Computation
nmax = 50
tol  = 1e-15

# Method to introduce initial conditions to refine (ignore if f1 = True)
# Same input used for f2 or f3 (if f2 is False)
method = 'insitu'       # Valid values: insitu / text

# If in situ, then enter manually on this section
x0  = 0.98883562973961558473
z0  = 0.00070003539306020212
vy0 = 0.00891271945118244681 # Recall y0 = vx0 = vz0 = 0 for a valid Halo
IC  = np.array([x0, z0, vy0])

# If text, then enter the .txt file on the current folder to retrieve the data from
# The favoured format is:
# It is accepted to have comments initiated by a # at any section of the text
# as long as the order of the values is consistent with the one given and no
# more than one set is given in the same .txt
# x0 = value # with no word after the values
# z0 = value
# vy0 = value
# although a full set of 6 ICs will be accepted and processed
# IC = 'sample.txt' # Uncomment this line to enable text

# Combining all input data into a single data structure
Input = {'flags': flags, 'tf': tf, 'nsteps': nsteps, 'opt': opt, 'LP': LP,
    'm': m, 'phi': phi, 'Az': Az, 'nmax': nmax, 'tol': tol, 'method': method,
    'IC': IC}

Halo_Main(Input)
