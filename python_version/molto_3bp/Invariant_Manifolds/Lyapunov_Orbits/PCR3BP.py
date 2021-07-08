def PCR3BP_propagator(S0, et0, deltat, prnt_out_dt, stop_fun, *params):
# ------------------- NUMERICAL INTEGRATION FUNCTION ----------------------
# This function propagates the state of a S/C according to the PBR3BP.
#
# -------------------------------------------------------------------------
#
# Import required functions
#
    from scipy.integrate import solve_ivp
    import numpy as np

# Initial S/C state vector in the input inertial reference frame
    state = S0
# Final ephemeris time
    etf = et0 + deltat

# Function that computes the state derivatives
    derivs = PCR3BP_state_derivs

# # Maximum integration time step is 1/50 of the orbital period
    max_step = 1e-3

# # Options of the Runge Kutta solver
    options = {'RelTol' : 5e-14, 'AbsTol' : 1e-30, 'MaxStep' : max_step}

# Avoid errors when the print out step has been set higher than the
# propagation duration
    if (prnt_out_dt > abs(etf-et0)):
        prnt_out_dt = etf-et0

# Set events function to the input function handle
    if stop_fun != 'None':
        options['Events'] = stop_fun
        try:
            if stop_fun.x:
                options['Events'].direction = -state[-1]
        except:
            pass
        options['Events'].terminal = True
# ---------- SOLVE FOR THE TRAJECTORY WITH AN ODE45 INTEGRATOR ------------
        sol = solve_ivp(derivs, (et0, etf),
            state, t_eval = np.linspace(et0, etf, int(abs((etf-et0)/prnt_out_dt)) +1),
            events = options['Events'],
            args = params,
            rtol = options['RelTol'],
            atol = options['AbsTol'],
            max_step = options['MaxStep'])

        times  = np.append(sol.t[:], sol.t_events[0][:], axis = 0)
        states = np.append(sol.y[:], sol.y_events[0][0][:, np.newaxis], axis = 1)
    else:
        sol = solve_ivp(derivs, (et0, etf),
            state, t_eval = np.linspace(et0, etf, int(abs((etf-et0)/prnt_out_dt)) +1),
            args = params,
            rtol = options['RelTol'],
            atol = options['AbsTol'],
            max_step = options['MaxStep'])

        times  = sol.t
        states = sol.y


# Update the final S/C state value and ephemeris time
    SF = states[:,-1]
    etf = times[-1]

    return SF, etf, states, times

def PCR3BP_state_derivs (t, state, mu1, mu2):
# This function computes the acceleration acting at time "et" on
# the fourth body according to the PBRFBP model.
# INPUTS:
#  - t     : Current non-dimensional time [sec]
#  - state  : Current non-dimensional S/C state in the synodic reference frame
#             [-,-]
# OUTPUTS:
#  - derivs: Total derivatives of the S/C state vector
#     - derivs (1:2) : Position derivatives
#     - derivs (3:4) : Velocity derivatives
##############################################################################
    import numpy as np

# Initialize derivative vector (4 components)
    derivs = np.zeros(4)
# Get state
    [x, y, xdot, ydot] = state
    # Get parameters (2D)
    Ux = mu2*(mu1 - x) - mu1*(mu2 + x)\
        - (mu2*(mu1 - x))/((mu1 - x)**2 + y**2)**(3/2)\
        + (mu1*(mu2 + x))/((mu2 + x)**2 + y**2)**(3/2)
    Uy = (mu1*y)/((mu2 + x)**2 + y**2)**(3/2) - mu2*y - mu1*y\
        + (mu2*y)/((mu1 - x)**2 + y**2)**(3/2)
    # Velocity derivatives (2D)
    derivs[2:] = [2*state[3] - Ux, -2*state[2] - Uy]

# Position derivatives
    derivs[:2] = state[2:]

    return derivs
