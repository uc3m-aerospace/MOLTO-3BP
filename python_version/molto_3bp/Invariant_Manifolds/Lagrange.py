def lagrange_points(mu, guess):
#
# Solve equation to obtain colinear Lagrange points
#
    from scipy.optimize import fsolve
    import numpy as np

    poslagrange1 = fsolve(findlagrange, guess[0], xtol = 1e-10, args = (mu,))
    poslagrange2 = fsolve(findlagrange, guess[1], xtol = 1e-10, args = (mu,))
    poslagrange3 = fsolve(findlagrange, guess[2], xtol = 1e-10, args = (mu,))

    pos = np.array([[poslagrange1[0], 0], [poslagrange2[0], 0], [poslagrange3[0], 0],
        [np.cos(np.pi/3)-mu, np.sin(np.pi/3)],
        [np.cos(np.pi/3)-mu, -np.sin(np.pi/3)]])

    return pos

def findlagrange(x, mustar):
    #This function contains the equation used to find the Lagrangian points.
    eqn = -((1-mustar)/abs(mustar+x)**3)*(mustar+x)+\
        (mustar/abs(1-mustar-x)**3)*(1-mustar-x) +x
    return eqn

def plot_lagrange_points(mu1, mu2, pos):

    import matplotlib.pyplot as plt
    plt.plot(-mu2, 0, 'ko')
    plt.plot(mu1, 0, 'ko')
    plt.plot(pos[0, 0], pos[0, 1], 'bo')
    plt.plot(pos[1, 0], pos[1, 1], 'bo')
    plt.plot(pos[2, 0], pos[2, 1], 'bo')
    plt.plot(pos[3, 0], pos[3, 1], 'ro')
    plt.plot(pos[4, 0], pos[4, 1], 'ro')
    plt.plot([-mu2, pos[3, 0]], [0, pos[3, 1]], '--')
    plt.plot([pos[3, 0], mu1], [pos[3, 1], 0], '--')
    plt.plot([-mu2, pos[4, 0]], [0, pos[4, 1]], '--')
    plt.plot([pos[4, 0], mu1], [pos[4, 1], 0], '--')

    moon_orbit = circle(-mu2, 0, 1)
    plt.plot(moon_orbit[:, 0], moon_orbit[:, 1], '--')
    plt.show()

def circle(x, y, r):
    # x and y: coordinates of the center of the circle
    # r: radius of the circle
    # 0.01 ~ angle step, bigger values will draw the circle faster but
    # imperfections may appear (not very smooth)
    import numpy as np

    ang  = np.linspace(0, 2*np.pi, int(2*np.pi/0.01 +1))
    xp   = r*np.cos(ang)
    yp   = r*np.sin(ang)
    circ = np.append(x+xp[:, np.newaxis], y+yp[:, np.newaxis], axis = 1)

    return circ
