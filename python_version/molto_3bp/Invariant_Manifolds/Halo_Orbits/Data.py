def Data(Data):
    Data['mE'] = 5.9722e24
    Data['mM'] = 7.342e22
    Data['mS'] = 1.98845114e30

    if Data['mode'] == 'SE':                      # Sun - Earth + Moon
        Data['mu'] = (Data['mE']+Data['mM'])/(Data['mE']+Data['mM']+Data['mS'])
    # Normalization relations
        Data['L']       = 1.497610041e6 # km
        Data['Ltheory'] = 1.496e8
        Data['TSE']     = 3.155815e7 # s
    elif Data['mode'] == 'EM':
        Data['mu']      = Data['mM']/(Data['mE']+Data['mM'])
        Data['L']       = 3.84388174e5
        Data['Ltheory'] = 3.84388174e5
        Data['TSE']     = 2.361e6 # s
    else:
        raise Exception('Halo Orbits:modeError.\
            The specified mode is outside range: [\'SE\', \'EM\']!')

    Data['params'] = (1-Data['mu'], Data['mu']) 

    return Data
