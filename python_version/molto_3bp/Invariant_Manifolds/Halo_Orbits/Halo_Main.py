def Halo_Main(Input):

    from .Data import Data

    Options = Data(Input)

    if (type(Options['flags'][0]) is not bool and type(Options['flags'][0]) is not int)\
        or (type(Options['flags'][1]) is not bool and type(Options['flags'][1]) is not int)\
        or (type(Options['flags'][2]) is not bool and type(Options['flags'][2]) is not int)\
        or (Options['flags'][0] == Options['flags'][1] and Options['flags'][0] == Options['flags'][2]\
        and Options['flags'][0] == 0):
        raise Exception('Halo_Main:FlagsError.\
            The values introduced for the flags are not valid!')

    if Options['flags'][0]:
        from .Halo_Generator import Halo_Generator
        if sum(Options['flags']) == 3:
            (Data, states_po, times_po, T_po, eigvec, eigval, inv_phi_0, xL) =\
                Halo_Generator(Options)
            return (Data, states_po, times_po, T_po, eigvec, eigval, inv_phi_0, xL)
        else:
            Halo_Generator(Options)
    else:
        if Options['method'] != 'insitu' and Options['method'] != 'text':
            raise Exception('Halo_Main:methodError.'+\
                '    The method selected is not valid [\'insitu\'][\'text\']!')
        if Options['flags'][1]:
            from .Halo_Num_Comp import Halo_Num_Comp
            Halo_Num_Comp(Options)
        elif Options['flags'][2]:
            from .Halo_Plot import Halo_Plot
            Halo_Plot(Options)
