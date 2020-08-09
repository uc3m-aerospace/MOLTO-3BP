###########################################################################
# LYAPUNOV TEST CASES
# The example consists on a Lyapunov orbit around L2 in the Earth - Moon System
# with a non-dimensional amplitude of .003 on the x axis
###########################################################################
#
# Import functions required
#
from Lyapunov_Orbits.Lyapunov import Lyapunov

# Main Functionalities (expected string)
# Accepted values: 'SE', 'EM'
f = 'EM'

# Case choice
Ax_tgt   = 3e-3 # Orbit's amplitude on the x axis

LP       = 2    # 1 -> Lagrange point L1
                # 2 -> Lagrange point L2

# Numerical parameters
prnt_out_dt = 0.01    # print time period

# Combining all input data into a single data structure
Input = {'mode': f, 'Ax_tgt': Ax_tgt, 'LP': LP, 'prnt_out_dt': prnt_out_dt}

Lyapunov(Input)
