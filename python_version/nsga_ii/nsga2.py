def nsga2(opt, *varargin):
# Function: result = nsga2(opt, varargin)
# Description: The main flowchart of of NSGA-II. Note:
#   All objectives must be minimization. If a objective is maximization, the
#   objective should be multipled by -1.
#
# Syntax:
#   result = nsga2(opt): 'opt' is generated by function nsgaopt().
#   result = nsga2(opt, param): 'param' can be any data type, it will be
#       pass to the objective function objfun().
#
#   Then ,the result structure can be pass to plotnsga to display the
#   population:  plotnsga(result);
#
# Parameters:
#   opt : A structure generated by funciton nsgaopt().
#   varargin : Additional parameter will be pass to the objective functions.
#       It can be any data type. For example, if you call: nsga2(opt, param),
#       then objfun would be called as objfun(x,param), in which, x is the
#       design variables vector.
# Return:
#   result : A structure contains optimization result.
#
# Import required functions
#
    import time
    from utilities.matrix_manipulation import value_matrix
    from nsga_ii.verifyOpt import verifyOpt
    from nsga_ii.evaluate import evaluate
    from nsga_ii.ndsort import ndsort
    from nsga_ii.statpop import statpop
    from nsga_ii.plotnsga import plotnsga
    from nsga_ii.callOutputfuns import callOutputfuns
    from nsga_ii.selectOp import selectOp
    from nsga_ii.crossover import crossover
    from nsga_ii.mutationOp import mutationOp
    from nsga_ii.extractPop import extractPop

#
# Start stopwatch to keep track of the execution times
#
    now = time.monotonic()

#*************************************************************************
# Verify the optimization model
#*************************************************************************
    opt = verifyOpt(opt)

#*************************************************************************
# variables initialization
#*************************************************************************
    nVar    = opt['numVar']
    nObj    = opt['numObj']
    nCons   = opt['numCons']
    popsize = opt['popsize']
    pop     = []

# pop : current population
# newpop : new population created by genetic algorithm operators
# combinepop = pop + newpop;
    for i in range(popsize):
        pop.append({'var' : value_matrix(1, nVar, 0.0),
            'obj' : value_matrix(1, nObj, 0.0),
            'cons' : value_matrix(1, nCons, 0.0),
            'rank' : 0,
            'distance' : 0,
            'prefDistance' : 0,     # preference distance used in R-NSGA-II
            'nViol' : 0,
            'violSum' : 0})

# state: optimization state of one generation
    state = {'current' : 1,     # current generation number
        'evaluateCount' : 0,    # number of objective function evaluation
        'totalTime' : 0,        # total time from the beginning
        'firstFrontCount' : 0,  # individual number of first front
        'frontCount' : 0,       # number of front
        'avgEvalTime' : 0}      # average evaluation time of objective function (current generation)

    result = {'opt' : opt, 'pops' : [], 'states' : []}
    for j in range(opt['maxGen']):
        result['pops'].append(pop)
        result['states'].append([state])

# global variables
    global STOP_NSGA
    STOP_NSGA = 0

#*************************************************************************
# initialize the P0 population
#*************************************************************************
    ngen = 1
    pop  = opt['initfun'][0](opt, pop, opt['initfun'][1:])
    [pop, state] = evaluate(opt, pop, state, varargin)
    [opt, pop]   = ndsort(opt, pop)

# state
    state['currentGen'] = ngen
    state['totalTime']  = time.monotonic()-now
    state               = statpop(pop, state)

    result['pops'][0]   = pop
    result['states'][0] = states

# output
    plotnsga(result, ngen)
    opt = callOutputfuns(opt, state, pop)

#*************************************************************************
# NSGA2 iteration
#*************************************************************************
    while ngen < opt['maxGen'] and STOP_NSGA == 0:
# 0. Display some information
        ngen += 1
        state['currentGen'] = ngen
        print('\n\n************************************************************')
        print(f"*      Current generation {ngen} / {opt['maxGen']}")
        print('************************************************************')

# 1. Create new population
        newpop = selectOp(opt, pop)
        newpop = crossover(opt, newpop, state)
        newpop = mutationOp(opt, newpop, state)
        [newpop, state] = evaluate(opt, newpop, state, varargin)

# 2. Combine the new population and old population : combinepop = pop + newpop
        combinepop = [pop, newpop]

# 3. Fast non dominated sort
        [opt, combinepop] = ndsort(opt, combinepop)

# 4. Extract the next population
        pop = extractPop(opt, combinepop)

# 5. Save current generation results
        state['totalTime'] = time.monotonic()-now
        state = statpop(pop, state)

        result['pops'][ngen-1] = pop
        result['states'][ngen-1]  = state

# 6. plot current population and output
        if ngen % opt['plotInterval'] == 0:
            plotnsga(result, ngen)
        if ngen % opt['outputInterval'] == 0:
            opt = callOutputfuns(opt, state, pop)

# call output function for closing file
        opt = callOutputfuns(opt, state, pop, -1)

        print('The complete running time is: ')
        print(time.monotonic()-now)
        return result
