def Data(Data):
    Data['mE'] = 5.9722e24
    Data['mM'] = 7.342e22
    Data['mS'] = 1.98845114e30

    if Data['opt'] == 1:                      # Sun - Earth + Moon
        Data['mu'] = (Data['mE']+Data['mM'])/(Data['mE']+Data['mM']+Data['mS'])
    # Normalization relations
        Data['L']       = 1.497610041e6 # km
        Data['Ltheory'] = 1.496e8
        Data['TSE']     = 3.155815e7 # s
    elif Data['opt'] == 2:
        Data['mu']      = Data['mM']/(Data['mE']+Data['mM'])
        Data['L']       = 3.84388174e5
        Data['Ltheory'] = 3.84388174e5
        Data['TSE']     = 2.361e6 # s
    else:
        raise Exception('Halo Orbits:OptError.\
            The specified value for the opt parameter is outside range: [1, 2]!')

    return Data
