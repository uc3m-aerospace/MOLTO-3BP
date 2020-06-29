##### HALO ORBITS GENERATION #####
#
# Importing required functions
#
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#### 1st approx - Linear equations solution

### Choice case
print("Input the option required by the user")
print("1 -> Sun - Earth+Moon System")
print("2 -> Earth - Moon Syst")
opt = int(input("opt = "))  # 1 -> Sun - Earth+Moon System
                            # 2 -> Earth - Moon Syst

print("Input the orientation of the orbit")
print("1 -> Nothern orbit   A_z > 0")
print("3 -> Southern orbit  A_z < 0")
m = int(input("m = "))      # 1 -> Nothern orbit   A_z > 0
                            # 3 -> Southern orbit  A_z < 0

phi = 0                 # rad  <-- Select initial value
Az  = 110e3             # km   <-- Select initial value
Ax  = 206e3             # km

t   = np.linspace(0, 100, int(1e4))   # s <-- Propagation time

# # DATA

mE = 5.9722e24
mM = 7.342e22
mS = 1.98845114e30

if opt==1:                      # Sun-Earth
    mu = (mE+mM)/(mE+mM+mS)     # Expected value 3.036e-6
elif opt == 2:
    mu = 1.215e-2
else:
    raise Exception('Halo Orbits:OptError.\
        The specified value for the opt parameter is outside range: [1, 2]!')

# Basic Variables Calculation

rh = (mu/3)**(1/3)
gamma1ch = rh*(1-1/3*rh-1/9*rh**2-1/27*rh**3)
gamma2ch = rh*(1+1/3*rh-1/9*rh**2)

# gamma polynomial
poly  = [1, -(3-mu), 3-2*mu, -mu, 2*mu, -mu] # for gamma 1 --> L1
poly2 = [1, (3-mu), 3-2*mu, -mu, -2*mu, -mu] # for gamma 2 --> L2

root1  = np.roots(poly)
root2  = np.roots(poly2)
gamma1 = np.real(root1[np.imag(root1) == 0])
gamma2 = np.real(root1[np.imag(root2) == 0])

# cn coefficients calculation
def cn1(n): return 1/gamma1**3*(+1**n*mu+(-1)**n*((1-mu)*gamma1**(n+1))/(1-gamma1)**(n+1)) # L1
def cn2(n): return 1/gamma2**3*(-1**n*mu+(-1)**n*((1-mu)*gamma2**(n+1))/(1+gamma2)**(n+1)) # L2

if opt == 1: # Check with Koon results
    print('---- CHECKs - with Koon ------\n')
    print('c_2 = %.10f' % cn1(2))
    print('c_3 = %.10f' % cn1(3))
    print('c_4 = %.10f\n' % cn1(4))

# Eigenvalues Calculation

c2 = cn1(2); c3 = cn1(3); c4 = cn1(4) # PROVISIONAL PARA COMPROBAR CON KOON

# lambda=sqrt((c2+sqrt(9*c2^2-8*c2))/2);
wp = np.sqrt((2-c2+np.sqrt(9*c2**2-8*c2))/2)
wv = np.sqrt(c2)
lmbda = wp

if opt == 1: # Check with Koon results
    print('w_p = %.10f' % wp)
    print('w_v = %.10f\n' % wv)

# Variables - Linear solution

kappa   = (wp**2+1+2*c2)/(2*wp)
kappa_C = 2*lmbda/(lmbda**2+1-c2)

# Third-Order Richardson Expansion - Intro

d1 = 3*lmbda**2/kappa*(kappa*(6*lmbda**2-1)-2*lmbda)

a21 = (3*c3*(kappa**2-2))/(4*(1+2*c2))
a22 = 3*c3/(4*(1+2*c2))
a23 = -3*c3*lmbda/(4*kappa*d1)*(3*kappa**3*lmbda-6*kappa*(kappa-lmbda)+4)
a24 = -3*c3*lmbda/(4*kappa*d1)*(2+3*kappa*lmbda)

b21 = -3*c3*lmbda/(2*d1)*(3*kappa*lmbda-4)
b22 = 3*c3*lmbda/d1

d21 = -c3/(2*lmbda**2)

s1 = (2*lmbda*(lmbda*(1+kappa**2)-2*kappa))**-1*\
    (3/2*c3*(2*a21*(kappa**2-2)-a23*(kappa**2+2)-2*kappa*b21)-3/8*c4*(3*kappa**4-8*kappa**2+8))
s2 = (2*lmbda*(lmbda*(1+kappa**2)-2*kappa))**-1*\
    (3/2*c3*(2*a22*(kappa**2-2)+a24*(kappa**2+2)+2*kappa*b22+5*d21)+3/8*c4*(12-kappa**2))


l1    = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-kappa**2)+2*lmbda**2*s1
l2    = 3/2*c3*(a24-2*a22)+9/8*c4+2*lmbda**2*s2
incre = wp**2-wv**2

if opt == 1: # Check with Koon results
    print('s_1 = %.10f' % s1)
    print('s_2 = %.10f' % s2)
    print('l_1 = %.10f' % l1)
    print('l_2 = %.10f' % s2)
    print('Increment = %.10f\n' % incre)

# Extra - for 3rd order Richardson

d2  = 8*lmbda**2/kappa*(kappa*(11*lmbda**2-1)-2*lmbda)

a31 = -9*lmbda/(4*d2)*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa**2))+\
    (9*lmbda**2+1-c2)/(2*d2)*(3*c3*(2*a23-kappa*b21)+c4*(2+3*kappa**2))
a32 = -1/d2*(9*lmbda/4*(4*c3*(kappa*a24-b22)+kappa*c4)+\
    3/2*(9*lmbda**2+1-c2)*(c3*(kappa*b22+d21-2*a24)-c4))

b31 = 3/(8*d2)*(8*lmbda*(3*c3*(kappa*b21-2*a23)-c4*(2+3*kappa**2))+\
    (9*lmbda**2+1+2*c2)*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa**2)))
b32 = 1/d2*(9*lmbda*(c3*(kappa*b22+d21-2*a24)-c4)+\
    3/8*(9*lmbda**2+1+2*c2)*(4*c3*(kappa*a24-b22)+kappa*c4))

d31 = 3/(64*lmbda**2)*(4*c3*a24+c4)
d32 = 3/(64*lmbda**2)*(4*c3*(a23-d21)+c4*(4+kappa**2))

# Extra

# Relationships

# Amplitude Constraint
AxC = np.sqrt(-(incre+l2*Az**2)/l1) # km
if AxC < Ax:
    print('Amplitude Constraint is Satisfied')
else:
    print('FAILED Amplitude Constraint!!')

# Phase Angle Relationship
if m in [1, 3]:
    psi = phi+m*np.pi/2 # rad
else:
    raise Exception('Halo Orbits:mError.\
        The specified value for the m parameter is outside range: [1 or 3]!')

# EoM

x = -Ax*np.cos(wp*t+phi)
y = kappa*Ax*np.sin(wp*t+phi)
z = Az*np.sin(wv*t+psi)

# # Plots

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.plot(x, y)
ax2.plot(x, z, 'tab:orange')
ax3.plot(y, z, 'tab:green')
for ax in (ax1, ax2, ax3):
    ax.set(xlabel='x', ylabel='y')
plt.show()
