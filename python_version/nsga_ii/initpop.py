def initpop(opt, pop, *varargin):
# Function: pop = initpop(opt, pop, varargin)
# Description: Initialize population.
# Syntax:
#   pop = initpop(opt, pop)
#     (default) Create a random initial population with a uniform distribution.
#
#   pop = initpop(opt, pop, 'pop.txt')
#     Load population from exist file and use the last population. If the popsize
#     less than the current popsize, then random numbers will used to fill the population.
#
#   pop = initpop(opt, pop, 'pop.txt', ngen)
#     Load population from file with specified generation.
#
#   pop = initpop(opt, pop, oldresult)
#     Specify exist result structure.
#
#   pop = initpop(opt, pop, oldresult, ngen)
#     Specify exist result structure and the population which will be used.
#
# Parameters:
#   pop : an empty population
#
#*************************************************************************
#*************************************************************************
# 1. Identify parameters
#*************************************************************************
    method = 'uniform'
    if (2 + len(varargin) >= 3):
        if varargin[0] is str:
            method = 'file'
        elif varargin[0] is dict:
            method = 'existpop'
#*************************************************************************
# 2. Initialize population with different methods
#*************************************************************************
    if method == 'uniform':
        pop = initpopUniform(opt,pop)
    elif method == 'file':
        print('...Initialize population from file ' + varargin[0] + '.')
        pop = initpopFromFile(opt, pop, varargin)
    elif method == 'existpop':
        print('...Initialize population from specified result.')
        pop = initpopFromExistResult(opt, pop, varargin)
    return pop

def initpopFromFile(opt, pop, *varargin):
# Function: pop = initpopFromFile(opt, pop, varargin)
# Description: Load population from specified population file.
# Syntax:
#   pop = initpop(opt, pop, 'pop.txt')
#   pop = initpop(opt, pop, 'pop.txt', ngen)
#
#*************************************************************************
#
# Import required functions
#
    from nsga_ii.loadpopfile import loadpopfile

    fileName = varargin[0]

    oldResult = loadpopfile(fileName)
    pop = initpopFromExistResult(opt, pop, oldResult, varargin[1:])
    return pop

def initpopFromExistResult(opt, pop, *varargin):
# Function: pop = initpopFromExistResult(opt, pop, varargin)
# Description: Load population from exist result structure.
# Syntax:
#   pop = initpop(opt, pop, oldresult)
#   pop = initpop(opt, pop, oldresult, ngen)
#
#*************************************************************************
#
# Import Warnings module for Python
#
    import warnings
# 1. Verify param
    oldresult = varargin[0]
    try:
        oldpops = oldresult['pops']
    except:
        raise Exception('NSGA2:InitPopError. The result structure specified is not correct!')

    if opt['numVar'] != len(oldpops[0][0]['var']) or\
        opt['numObj'] != len(oldpops[0][0]['obj']) or\
        opt['numCons'] != len(oldpops[0][0]['cons']):
        raise Exception('NSGA2:InitPopError.\
            The specified optimization result is not for current optimization model!')

# 2. Determine which population would be used
    ngen = 0
    if (2 + len(varargin) >= 4):
        ngen = varargin[1]

    maxGen = len(varargin)

    if ngen == 0:
        ngen = maxGen
    elif ngen > maxGen:
        warnings.warn(f'NSGA2:InitPopWarning.\
            The specified generation "{ngen}" does not exist, using "{maxGen}" instead.')

    ngen = maxGen

# 3. Create initial population
    popsizeOld = len(oldpops[0])
    popsizeNew = opt['popsize']

    if popsizeNew <= popsizeOld:    # a) All from old pop
        for i in range(0, popsizeNew):
            pop[i]['var'] = oldpops[ngen][i]['var']
    else:                           # b) Use random individuals to fill the population
        for i in range(0, popsizeOld):
            pop[i]['var'] = oldpops[ngen][i]['var']
        pop[popsizeOld:] = initpopUniform(opt, pop[popsizeOld:])
    return pop

def initpopUniform(opt, pop):
# Function: pop = initpopUniform(opt, pop)
# Description: Initialize population using random number
#
#*************************************************************************
#
# Import required functions
#
    from utilities.matrix_manipulation import rand_value_matrix, transpose, vector_multiply, matrix_subtraction
    from nsga_ii.varlimit import varlimit

    nVar = opt['numVar']
    type = opt['vartype']

    lb = opt['lb']
    ub = opt['ub']

    popsize = len(pop)
    for i in range(0, popsize):
        var1 = matrix_subtraction(ub, transpose(lb))
        var = vector_multiply(rand_value_matrix(1, nVar), var1[0])
        for j in range(len(var)):
            for k in range(len(var[0])):
                var[j][k] += lb[j][0]


# if desing variable is integer, round to the nearest integer
    print(var)
    for j in range(len(var)):
        for k in range(len(var[0])):
            if( type[j] == 2):
                var[j][k] = round(var[j][k])

# limit in the lower and upper bound
    var = varlimit(var, lb, ub)

    pop[i]['var'] = var
    return pop
