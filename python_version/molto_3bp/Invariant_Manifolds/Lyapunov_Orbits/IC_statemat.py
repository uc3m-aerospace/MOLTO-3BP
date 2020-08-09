def IC_statemat(Ax, a, b):

    import numpy as np
    from scipy import linalg

    # Compute state transition matrix
    A = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [a, 0, 0, 2], [0, -b, -2, 0]])

    [eigval, eigvec] = linalg.eig(A)
    Mat = eigvec[0, 2] + eigvec[0, 3]
    x0 = Ax # Initial position in the x axis at a distance x0 from Libration point
    y0 = 0

    vx0 = x0/Mat*(eigvec[2,2] + eigvec[2,3])
    vy0 = x0/Mat*(eigvec[3,2] + eigvec[3,3])

    inv_phi_0 = linalg.inv(eigvec)

    return [x0, y0, vx0, vy0, eigvec, eigval, inv_phi_0]
