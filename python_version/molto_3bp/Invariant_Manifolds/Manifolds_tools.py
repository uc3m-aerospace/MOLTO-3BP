def construct(params, T0, states_po, times_po, eigvec, eigval, inv_phi_0,
    prnt_out_dt, npoints, stop_fun):

    #
    # Import required functions
    #
    import numpy as np
    from scipy import linalg
    from PCR3BP import PCR3BP_propagator

    ## Step 4: Computation of invariant manifolds
    T_po = T0
    idx  = np.linspace(0, len(times_po), npoints, endpoint = False, dtype = int)
    et0  = 0
    eps  = 1e-7
    deltat = 5*T_po
    sign   = np.array([-1, 1])

    # Matrix initialization
    states_u = []
    SF_u     = []
    times_u  = []
    states_s = []
    times_s  = []
    SF_s     = []

    for i in range(npoints):
        for j in range(len(sign)):
            x = states_po[:, idx[i]]
            t = times_po[idx[i]]
            phi = [eigvec[:, 0]*np.exp(eigval[0]*t), eigvec[:, 1]*np.exp(eigval[1]*t),
                eigvec[:, 2]*np.exp(eigval[2]*t), eigvec[:, 3]*np.exp(eigval[3]*t)]
            Phi = phi * inv_phi_0
            [mon_eigvals, mon_eigvecs] = linalg.eig(Phi)

            Yu = mon_eigvecs[:, 0] # unstable eigenvector
            Ys = mon_eigvecs[:, 1] # stable eigenvector

            Xu = np.real(states_po[:, idx[i]] + eps*sign[j]*Yu)
            Xs = np.real(states_po[:, idx[i]] + eps*sign[j]*Ys)

            # Integrate unstable manifold forwards in time
            [SF, etf_u, states, times] = PCR3BP_propagator(Xu, params, et0, deltat,
                prnt_out_dt, stop_fun)

            states_u.append(states)
            times_u.append(times)
            SF_u.append(SF)

            # Integrate stable manifold backwards in time
            [SF, etf, states, times] =  PCR3BP_propagator(Xs, params, et0, -deltat,
                prnt_out_dt, stop_fun)

            states_s.append(states)
            times_s.append(times)
            SF_s.append(SF)

            print('Iteration: ' + str(2*i + j))

    return [states_s, times_s, SF_s, states_u, times_u, SF_u]

def plotm(mu1, mu2, pos, states_po, states_s, SF_s, states_u):

    import matplotlib.pyplot as plt
    import numpy as np

    plt.plot(mu1, 0, 'ko')
    plt.plot(pos[1, 0], pos[1, 1], 'ko')
    plt.plot(states_po[0], states_po[1], 'k')
    plt.show()

    for i in range(len(states_u)):
        plt.plot(states_u[i][0], states_u[i][1], 'r')
        plt.plot(states_s[i][0], states_s[i][1], 'b')
    # plt.plot(-mu2, 0 , 'bo')
    plt.plot(mu1, 0 , 'ro')
    plt.plot(pos[1, 0], pos[1, 1], 'ko')
    plt.plot(states_po[0], states_po[1], 'k')
    plt.plot(np.array([mu1-2.5e-3, mu1]), np.zeros(2), 'g--')
    plt.show()

    # plt.plot(-mu2, 0 , 'bo')
    plt.plot(mu1, 0 , 'ro')
    plt.plot(pos[:2, 0], pos[:2, 1], 'ko')
    plt.show()

    for j in SF_s:
        plt.plot(j[1], j[3], 'o')
    plt.show()
