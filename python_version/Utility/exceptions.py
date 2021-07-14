from Common.Service.benchmark import Benchmark as b
from Common.Config.configurator import Configurator as c


class UnknownKernelType(Exception):
    """
    Exception for unknown Kernel types
    """
    def __init__(self, message='Could not find kernel file gm_de431.tpc!!', kernel_file='gm_de431.tpc'):
        super().__init__(message)
        b.stat(f'Could not load kernel file {kernel_file}')


## Halo Orbits Exceptions ##

class HaloNumCompError(Exception):
    """
    Exception for Halo Numerical Computation Error
    """
    def __init__(self):
        super().__init__('Halo_Num_Comp:ICError.\n' + \
        'The text file selected does not have the right format!')

class HaloOrbitLPError(Exception):
    """
    Exception for Halo Numerical Computation Error
    """
    def __init__(self):
        super().__init__('Halo Orbits:LPError.\n\
                The specified value for the LP parameter is outside range: [1, 2]!')

class HaloOrbitmError(Exception):
    """
    Exception for Halo Numerical Computation Error
    """
    def __init__(self):
        super().__init__('Halo Orbits:mError.\
                The specified value for the m parameter is outside range: [1 or 3)!')

class HaloMainFlagsError(Exception):
    """
    Exception for Halo Numerical Computation Error
    """
    def __init__(self):
        super().__init__('Halo_Main:FlagsError.\
                The values introduced for the flags are not valid!')

class HaloMainMethodError(Exception):
    """
    Exception for Halo Numerical Computation Error
    """
    def __init__(self):
        super().__init__('Halo_Main:methodError.' + \
                                '    The method selected is not valid [\'insitu\'][\'text\']!')




## ----------
## Lypaunov Orbits Exceptions ##


## ---------
## Manifold Sequence/Tools Exceptions

class ManifoldsSequenceError(Exception):
    """
    Exception for Manifolds Sequence Error
    """

    def __init__(self):
        super().__init__('Manifolds_Seq:typeError.'+\
                'The type selected is not valid [\'LY\'][\'HL\']!')

