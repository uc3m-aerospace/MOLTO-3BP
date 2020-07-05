##### HALO ORBITS PLOTTING TOOL #####
#
# Importing required functions
#
import numpy as np
from scipy.integrate import solve_ivp
from intFun import ThreeBodyProp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Use Results from Main2_Halo_num.m
# Initial conditions

x0  = 1.1107439873579033
y0  = 0
z0  = 0.035680331960522345
vx0 = 0
vy0 = 0.20363565695006591
vz0 = 0

q0 = np.array([x0, y0, z0, vx0, vy0, vz0])

tf = 3.39 # Period Time Insert
tspan = np.linspace(0, tf, int(1e5))
sol = solve_ivp(ThreeBodyProp, [0, tf],
    q0, t_eval = tspan)
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
