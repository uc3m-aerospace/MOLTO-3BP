def load_variables():
    #
    # Import required functions
    #
    import spiceypy as spice
    from kernels.Load_Kernels import load_kernels

    ## Load Spice kernels
    load_kernels()

    ## INITIAL DATA

    # Gravitational constants:
    mu_earth = spice.bodvrd('EARTH', 'GM', 1)[1][0] # Earth grav. pot [km3/s2]
    mu_sun   = spice.bodvrd('SUN', 'GM', 1)[1][0]   # Sun grav. pot [km3/s2]
    mu_moon  = spice.bodvrd('MOON', 'GM', 1)[1][0]  # Moon grav. pot [km3/s2]

    # Earth-Moon system
    mu_EM  = mu_moon/(mu_moon+mu_earth)

    L_EM = 3.85e5  # distance between primaries (p 33)
    T_EM = 2.361e6 # orbital period of primaries 27.32 days

    # Sun-Earth system
    mu_SE  = (mu_earth+mu_moon)/(mu_earth+mu_sun+mu_moon)

    L_SE = 1.496e8 # Distance Sun-Earth
    T_SE = 3.156e7 # Period of Earth around Sun: 365 days

    params_EM = {'mu1': (1-mu_EM), 'mu2': mu_EM, 'L_EM': L_EM, 'T_EM': T_EM}
    params_SE = {'mu1': (1-mu_SE), 'mu2': mu_SE, 'L_SE': L_SE, 'T_SE': T_SE}

    return [params_EM, params_SE]
