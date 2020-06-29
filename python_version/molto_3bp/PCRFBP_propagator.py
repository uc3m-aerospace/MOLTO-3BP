def PCRFBP_propagator(S0, params, et0, deltat, prnt_out_dt, stop_fun):
# ------------------- NUMERICAL INTEGRATION FUNCTION ----------------------
# This function propagates the state of a S/C according to the PBRFBP.
#
# -------------------------------------------------------------------------
#
# Import required functions
#
    from molto_3bp.PCRFBP_state_derivs import PCRFBP_state_derivs
    from scipy.integrate import solve_ivp
    import numpy as np
# Initial S/C state vector in the input inertial reference frame
    state = S0
# Final ephemeris time
    etf = et0 + deltat

# Function that computes the state derivatives
    derivs = PCRFBP_state_derivs

# # Options of the Runge Kutta solver
# # Maximum integration time step is 1/50 of the orbital period
    max_step = 1e20 # 1/50 Moon orbit of 100 km altitude

# Additional exit condition of the Cowell's propagator (if it is reached
# before the final propagation time)
    options = {'RelTol' : 1e-9, 'AbsTol' : 1e-9, 'MaxStep' : max_step}

# Avoid errors when the print out step has been set higher than the
# propagation duration
    if (prnt_out_dt > etf-et0):
        prnt_out_dt = etf-et0

# Set events function to the input function handle
    if stop_fun != 'None':
        options['Events'] = stop_fun

# ---------- SOLVE FOR THE TRAJECTORY WITH AN ODE45 INTEGRATOR ------------
        sol = solve_ivp(derivs, np.linspace(et0, etf, int((etf-et0)/prnt_out_dt)),
            state, events = options['Events'],
            rtol = options['RelTol'],
            atol = options['AbsTol'],
            max_step = options['MaxStep'])
    else:
        sol = solve_ivp(derivs, np.linspace(et0, etf, int((etf-et0)/prnt_out_dt)),
            state,
            rtol = options['RelTol'],
            atol = options['AbsTol'],
            max_step = options['MaxStep'])

# When propagating in time, assure that the final propagation time is
# included in the simulation
    if (stop_fun == 'None' and times[-1] < etf):
        time0 = times[-1]
        timef = etf
        S0 = states[-1]
        sol = solve_ivp(derivs, np.linspace(time0, timef, 2), S0,
            rtol = options['RelTol'],
            atol = options['AbsTol'],
            max_step = options['MaxStep'])
        times = np.append(times, ttt[-1])
        states = np.append(states, sss[:,-1,np.newaxis], axis = 1)

# Update the final S/C state value and ephemeris time
    SF = states[:,-1]
    etf = times[-1]

    return SF, etf, states, times
