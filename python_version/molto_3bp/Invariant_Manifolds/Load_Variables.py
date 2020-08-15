def load_variables(Data):
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

    if Data['mode'] == 'SE':                        # Earth-Moon system
        Data['mu'] = (mu_earth+mu_moon)/(mu_earth+mu_sun+mu_moon)
    # Normalization relations
        Data['L']       = 1.497610041e6 # km        # Distance Sun-Earth
        Data['Ltheory'] = 1.496e8
        Data['TSE']     = 3.155815e7 # s            # Period of Earth around Sun: 365 days
    elif Data['mode'] == 'EM':                      # Sun-Earth system
        Data['mu']      = mu_moon/(mu_moon+mu_earth)
        Data['L']       = 3.84388174e5              # distance between primaries
        Data['Ltheory'] = 3.84388174e5
        Data['TSE']     = 2.361e6 # s               # orbital period of primaries 27.32 days
    else:
        raise Exception('Halo Orbits:modeError.\
            The specified mode is outside range: [\'SE\', \'EM\']!')

    Data['params'] = (1-Data['mu'], Data['mu'])

    return Data
