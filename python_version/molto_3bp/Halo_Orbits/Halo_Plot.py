def Halo_Plot(Data):
    ##### HALO ORBITS PLOTTING TOOL #####
    #
    # Importing required functions
    #
    import numpy as np
    from scipy.integrate import solve_ivp
    from intFun import ThreeBodyProp
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Initial conditions
    if Data['flags'][0] or Data['flags'][1] or Data['method'] == 'insitu':
        [x0, z0, vy0] = Data['IC']
    else:
        fid  = open(Data['IC'],'r')
        info = fid.readlines()
        IC = []
        for i in info:
            if i.split()[0] == '#':
                IC.append(i.split()[-1])
        if len(IC) == 3:
            [x0, z0, vy0] = IC
        elif len(IC) == 6:
            [x0, z0, vy0] = [IC[0], IC[2], IC[4]]
        else:
            raise Exception('Halo_Num_Comp:ICError.' +\
                '   The text file selected does not have the right format!')

    q0 = np.array([x0, 0, z0, 0, vy0, 0])
    tspan = np.linspace(0, Data['tf'], int(Data['nsteps']))

    sol = solve_ivp(ThreeBodyProp, [0, Data['tf']],
        q0, t_eval = tspan, args = (Data['mu'],),
        atol = 1e-10,
        rtol = 1e-7)
    t = sol.t
    q = sol.y
    x = q[0]; y = q[1]; z = q[2]

    # Figures

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(x, y, z)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(x, y)
    ax1.set(xlabel='x', ylabel='y')
    ax2.plot(x, z)
    ax2.set(xlabel='x', ylabel='z')
    ax3.plot(y, z)
    ax3.set(xlabel='y', ylabel='z')
    plt.show()
