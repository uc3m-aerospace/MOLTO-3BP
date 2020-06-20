def varlimit(var, lb, ub):
# Function: var = varlimit(var, lb, ub)
# Description: Limit the variables in [lb, ub].
#
#*************************************************************************
    numVar = len(var)
    for i in range(0, numVar):
        if var[i][0] < lb[i][0]:
            var[i][0] = lb[i][0]
        elif var[i][0] > ub[0][i]:
            var[i][0] = ub[0][i]
