#--------------------------------------------------------------------------
def molto_3bp(input):
#--------------------------------------------------------------------------
#	MOLTO-3BP Software Computation Core
#
#	This program is developed at the Universidad Carlos III de Madrid.
#
#   The program is released under the MIT License
#
#--------------------------------------------------------------------------
#
#    Function that defined the parameters for NSGA-II and call the
#    required routines to solve the interplanetary transfer
#
#--------------------------------------------------------------------------
#
# Import mathematical packages for consts and operations
#
    import math
    from utilities.matrix_manipulation import value_matrix
#
# Import additional functions required
#
    from nsga_ii.nsga2 import nsga2
    from nsga_ii.nsgaopt import nsgaopt
    from molto_3bp.fitness_nsga2 import fitness_nsga2
    from nsga_ii.initpop import initpop
#
# Define constants
#
# Gravitational constants:
#
    mu_sun   = 1.32712440018e11   # Sub grav. pot [km3/s2]
    mu_earth = 3.986004418e5      # Earth grav. pot [km3/s2]
    mu_moon  = 4.9048695e3        # Moon grav. pot [km3/s2]
#
# Non-dimensional constants (Simó et al. (1995)):
#
    mu    = mu_moon/(mu_earth+mu_moon)      # Earth-Moon mass parameter [-]
    m_s   = mu_sun/(mu_earth+mu_moon)       # Scaled mass of the Sun
    rho   = 3.88811143*10**2                # Scaled Sun-(Earth+Moon) distance
    om_s  = -9.25195985*10**(-1)            # Scaled angular velocity of the Sun
    l_em  = 3.84405000*10**8                # [m] Earth-Moon distance
#
    input['rho']  = rho
    input['mu']   = mu
    input['m_s']  = m_s
    input['om_s'] = om_s
#
    input['DU']    = l_em                   # [m] Distance unit
    input['TU']    = 4.34811305*24*3600     # [days] Time unit
    input['VU']    = 1.02323281*10**3       # [m/s] Speed unit
#
    input['r0'] = input['r0']*10**3/input['DU'];
    input['rf'] = input['rf']*10**3/input['DU'];
#
#
# Set Min/Max Values for the population
#
    LB = [[0.0], [1.0], [0.0], [1.0]]
    UB = [[2*math.pi, math.sqrt(2), 2*math.pi/(om_s), 30.0]]
#
# Set the vartype (=1 Continuos) (=2 Discrete)
#
    vartype = value_matrix(1,len(LB),1.0)
#
# Set the number of variables
#
    nvars = len(LB)
#
#--------------------------------------------------------------------------
# Set the NSGA-II genetic algorithm parameters
#--------------------------------------------------------------------------
#
#  Initialize options structure
#
    options = nsgaopt()
#
# Set User defined population size
#
    if 'popsize' in input:
        options['popsize'] = input['popsize']
#
# Set User defined Maximum number of generations
#
    if 'maxGen' in input:
        options['maxGen'] = input['maxGen']
#
    options['numObj']  = 2       # number of objectives
    options['numVar']  = nvars   # number of design variables
    options['numCons'] = 0       # number of constraints
    options['lb']      = LB      # lower bound of x
    options['ub']      = UB      # upper bound of x
#
    options['nameObj']    = ['Delta r (km)','Delta V (km/s)']  # List of objectives names for the plot
    options['crossover']  = ['intermediate', 0.5]   # crossover operator -- Intermediate crossover [3]
# creates two children from two parents: parent1 and parent2. If it lies in the range [0, 1],
# the children created are within the two parent. If algorithm is premature, try to set ratio larger than 1.0.
    options['mutation']   = ['gaussian',0.1, 0.2]   # mutation operator
# (scale = 0.1 (deviation of the random number), shrink=0.5)
# --> for example, shrink inside [0.5, 1.0] is usually used for local search.
# A large mutation range (shrink == 0) is require getting out of the local Pareto-optimal fronts
    options['crossoverFraction'] = 0.8  # crossover fraction of variables of an individual ( 2/numVar )
# --> only crossoverFraction of all variables would do crossover
    options['mutationFraction']  = 0.3  # only mutaionFraction of all variables would do mutation (default 2/numvar)
    options['objfun']            = fitness_nsga2    # objective function handle
    options['plotInterval']      = 1    # interval between two calls of "plotnsga".
    options['outputInterval']    = 1
    options['outputfile']        = input['output_file']    # outputfile
    options['useParallel']       = input['useParallel']    # Parallel option.
    options['vartype']           = vartype
#
    if (len(input['init_file'])!=0):
        options['initfun'] =  [initpop, input['init_file']]
#
#--------------------------------------------------------------------------
# Call NSGA-II algorithm
#--------------------------------------------------------------------------
#
    output = nsga2(options,input)
#
# Direct results to input (Main)
#
    return output
