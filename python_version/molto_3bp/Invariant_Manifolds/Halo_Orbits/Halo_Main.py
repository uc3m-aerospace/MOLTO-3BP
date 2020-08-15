def Halo_Main(Input):

    if (type(Input['flags'][0]) is not bool and type(Input['flags'][0]) is not int)\
        or (type(Input['flags'][1]) is not bool and type(Input['flags'][1]) is not int)\
        or (type(Input['flags'][2]) is not bool and type(Input['flags'][2]) is not int)\
        or (Input['flags'][0] == Input['flags'][1] and Input['flags'][0] == Input['flags'][2]\
        and Input['flags'][0] == 0):
        raise Exception('Halo_Main:FlagsError.\
            The values introduced for the flags are not valid!')

    if Input['flags'][0]:
        from .Halo_Generator import Halo_Generator
        if sum(Input['flags']) == 3:
            (states_po, times_po, T_po, eigvec, eigval, inv_phi_0, xL) =\
                Halo_Generator(Input)
            return (states_po, times_po, T_po, eigvec, eigval, inv_phi_0, xL)
        else:
            Halo_Generator(Input)
    else:
        if Input['method'] != 'insitu' and Input['method'] != 'text':
            raise Exception('Halo_Main:methodError.'+\
                '    The method selected is not valid [\'insitu\'][\'text\']!')
        if Input['flags'][1]:
            from .Halo_Num_Comp import Halo_Num_Comp
            Halo_Num_Comp(Input)
        elif Input['flags'][2]:
            from .Halo_Plot import Halo_Plot
            Halo_Plot(Input)
