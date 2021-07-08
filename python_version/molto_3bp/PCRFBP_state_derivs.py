def PCRFBP_state_derivs (t, state, params):
# This function computes the acceleration acting at time "et" on
# the fourth body according to the PBRFBP model.
# INPUTS:
#  - et     : Current non-dimensional time [sec]
#  - state  : Current non-dimensional S/C state in the synodic reference frame
#             [-,-]
# OUTPUTS:
#  - derivs: Total derivatives of the S/C state vector
#     - derivs (1:2) : Position derivatives
#     - derivs (3:4) : Velocity derivatives
##############################################################################
# Initialize derivative vector (4 components)
    derivs = np.zeros(4)
# Get state
    x    = state[0]
    y    = state[1]
    xdot = state[2]
    ydot = state[3]

# Get parameters
    rho  = params['rho']
    mu   = params['mu']
    m_s  = params['m_s']
    om_s = params['om_s']

    dOm4dx = x - (m_s*np.cos(om_s*t))/rho**2 - (mu*(mu + x - 1))/((mu + x - 1)**2 + y^2)**(3/2)\
        + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)**2 + y**2)**(3/2))\
        - (m_s*(2*x - 2*rho*np.cos(om_s*t)))/(2*((x - rho*np.cos(om_s*t))**2\
        + (y - rho*np.sin(om_s*t))**2)**(3/2))
    dOm4dy = y - (m_s*np.sin(om_s*t))/rho**2 - (mu*y)/((mu + x - 1)**2 + y**2)**(3/2)\
        - (m_s*(2*y - 2*rho*np.sin(om_s*t)))/(2*((x - rho*np.cos(om_s*t))**2\
        + (y - rho*np.sin(om_s*t))**2)**(3/2)) + (y*(mu - 1))/((mu + x)**2 + y**2)**(3/2)

# Position derivatives
    derivs[:1] = state[2:]

# Velocity derivatives
    derivs[2:] = [dOm4dx + 2*state(4), dOm4dy - 2*state(3)])

    return derivs
