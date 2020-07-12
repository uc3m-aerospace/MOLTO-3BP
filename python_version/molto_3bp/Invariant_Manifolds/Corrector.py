def Corrector(Fun, x0, params, t0, Itmax, Tol, TolRel, TolAbs, dh, Ind_Fix):
# Periodic Orbit Corrector for autonomous system
# Fun                       - > function F = Fun(t,x) that gives the
#                                right-hand side
# x0                        - > Initial Condition guess
# t0                        - > Period guess
# Itmax, Tol, TolRel,TolAbs - > Numerical Parameters
# dh                        - > Jacobian: Finite differences step
# Ind_Fix                   - > (Integer) Index of the state vector variable
#                                that is fixed

###########################################################################
    #
    # Import required functions
    #
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy import linalg

    N    = len(x0)
    I    = np.identity(N)
    xFix = x0[Ind_Fix]
    for i in range(1, Itmax+1):

        # Compute the Initial Error
        x0[Ind_Fix] = xFix
        sol = solve_ivp(Fun, (0, t0), x0,
            args = params,
            rtol = TolRel,
            atol = TolAbs)
        t = sol.t; x = sol.y
        error = max(abs(x[:, 0] - x[:, -1]))
        print('Corrector: iteration = ' + str(i) + '  Error = %.5f' % error)

        xvar = np.append(x0, I.reshape(1, -1))

        if error < Tol:
            if i == 1:
                sol = solve_ivp(variational, (0, t0), xvar,
                    args = (params[0], params[1], Fun, N, dh),
                    rtol = TolRel,
                    atol = TolAbs)
                t = sol.t; x = sol.y
                M = x[N:, -1].reshape(N, -1, order = 'F')
                [Floquet, vec] = linalg.eig(M)
            break

        # Compute Monodromy Matrix
        sol = solve_ivp(variational, (0, t0), xvar,
            args = (params[0], params[1], Fun, N, dh),
            rtol = TolRel,
            atol = TolAbs)
        t = sol.t; x = sol.y
        xf = x[:, -1]
        M = x[N:, -1].reshape(N, -1, order = 'F')

        # Compute the derivative
        df = Fun(t0, x0, params[0], params[1])

        # Prepare equation A * [DeltT X1 X2... XN ] = b (without X(Index))
        A             = M - I
        A[:, Ind_Fix] = df
        B             = -(xf[:N] - x0)
        correc        = linalg.solve(A, B)
        for j in range(N):
            if j == Ind_Fix:
                t0 = t0 + correc[j]
            else:
                x0[j] = x0[j] + correc[j]
        [Floquet, vec] = linalg.eig(M)
        if i == Itmax:
            sol = solve_ivp(Fun, (0, t0), x0,
                args = params,
                rtol = TolRel,
                atol = TolAbs)
            t = sol.t; x = sol.y
            error = max(abs(x[:, 0] - x[:, -1]))
            print('Corrector: iteration = ' + str(i) + '  Error = %.5f' % error)

def variational(t, YV, mu1, mu2, Fun, N, dh):

    import numpy as np
    Jac = Jac_num(Fun, t, YV[:N], dh, mu1, mu2)
    YVtemp  = YV[N:].reshape(N, -1, order = 'F')
    deltadf = (Jac @ YVtemp).reshape(1, -1)[0]
    df  = np.append(Fun(t, YV[:N], mu1, mu2), deltadf)

    return df

def Jac_num(function, t, Y, dh, mu1, mu2):

    import numpy as np

    a = len(Y)
    Jnum = np.zeros((a, a))

    for i in range(a):
        YM = Y
        Ym = Y

        YM[i] = YM[i] + dh/2
        fM    = function(t, YM, mu1, mu2)

        Ym[i] = Ym[i] - dh/2
        fm    = function(t, Ym, mu1, mu2)

        Jnum[:, i] = (fM - fm)/dh

    return Jnum
