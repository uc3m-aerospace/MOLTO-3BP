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

print("Input the Lagrange Point analyzed")
print("1 -> Lagrange point L1")
print("2 -> Lagrange point L2")
LP = int(input("LP = "))  # 1 -> Sun - Earth+Moon System
                            # 2 -> Earth - Moon Syst

print("Input the orientation of the orbit")
print("1 -> Nothern orbit   A_z > 0")
print("3 -> Southern orbit  A_z < 0")
m = int(input("m = "))      # 1 -> Nothern orbit   A_z > 0
                            # 3 -> Southern orbit  A_z < 0

phi = 0                     # rad  <-- Select initial value
AzDim  = 220e3              # km   <-- Select initial value
AxDim  = 206e3              # km

t   = np.linspace(0, 100, int(1e4))   # s <-- Propagation time

# # DATA

mE = 5.9722e24
mM = 7.342e22
mS = 1.98845114e30

if opt==1:                      # Sun-Earth
    mu = (mE+mM)/(mE+mM+mS)     # Expected value 3.036e-6
# Normalization relations
# Intro to Lpoints
    L       = 1.497610041e6 # 1.845287088696577e+06; # km 1.428571428571428e+06 # REVISAR!!!!!!!!!!!!!!
    Ltheory = 1.496e8
    TSE     = 3.155815e7 # 3.2035e7; # s # REVISAR!!!!!!!!!!!!!!
elif opt == 2:
    mu      = 0.012150585609624 # 1.215e-2;
    L       = 3.84388174e5 # REVISAR!!!!!!!!!!!!!!
    Ltheory = 3.84388174e5
    TSE     = 2.361e6 # s # REVISAR!!!!!!!!!!!!!!
else:
    raise Exception('Halo Orbits:OptError.\
        The specified value for the opt parameter is outside range: [1, 2]!')
Tconversion = TSE/(2*np.pi)

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
gamma2 = np.real(root2[np.imag(root2) == 0])

# cn coefficients calculation
cn1 = lambda n: 1/gamma1**3*(+1**n*mu+(-1)**n*((1-mu)*gamma1**(n+1))/(1-gamma1)**(n+1)) # L1
cn2 = lambda n: 1/gamma2**3*(-1**n*mu+(-1)**n*((1-mu)*gamma2**(n+1))/(1+gamma2)**(n+1)) # L2

if LP == 1:
    c2 = cn1(2); c3 = cn1(3); c4 = cn1(4) # PROVISIONAL PARA COMPROBAR CON KOON
elif LP == 2:
    c2 = cn2(2); c3 = cn2(3); c4 = cn2(4) # PROVISIONAL PARA COMPROBAR CON KOON
else:
    raise Exception('Halo Orbits:LPError.\
        The specified value for the LP parameter is outside range: [1, 2]!')

r1 = Ltheory*((1-mu)-gamma1) # REVISAR!!!!!!!!!!!!
r2 = Ltheory*((1-mu)+gamma2)

if opt == 1: # Check with Koon results
    print('---- CHECKs - with Koon ------\n')
    print('c_2 = %.10f' % c2)
    print('c_3 = %.10f' % c3)
    print('c_4 = %.10f\n' % c4)

# Eigenvalues Calculation

# lambda = sqrt((c2+sqrt(9*c2^2-8*c2))/2)
wp = np.sqrt((2-c2+np.sqrt(9*c2**2-8*c2))/2)
wv = np.sqrt(c2)
lmbda = wp # For halo orbits

if opt == 1: # Check with Koon results
    print('w_p = %.10f' % wp)
    print('w_v = %.10f\n' % wv)

# Variables - Linear solution

kappa   = (wp**2+1+2*c2)/(2*wp)
kappa_C = 2*lmbda/(lmbda**2+1-c2)

# Third-Order Richardson Expansion

d1 = 3*lmbda**2/kappa*(kappa*(6*lmbda**2-1)-2*lmbda)
d2 = 8*lmbda**2/kappa*(kappa*(11*lmbda**2-1)-2*lmbda)

a21 = (3*c3*(kappa**2-2))/(4*(1+2*c2))
a22 = 3*c3/(4*(1+2*c2))
a23 = -3*c3*lmbda/(4*kappa*d1)*(3*kappa**3*lmbda-6*kappa*(kappa-lmbda)+4)
a24 = -3*c3*lmbda/(4*kappa*d1)*(2+3*kappa*lmbda)

b21 = -3*c3*lmbda/(2*d1)*(3*kappa*lmbda-4)
b22 = 3*c3*lmbda/d1

d21 = -c3/(2*lmbda**2)

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

# Relationships

# Amplitude Constraint

Ax = AxDim/L; Az = AzDim/L
AxC = np.sqrt(-(incre+l2*Az**2)/l1) # km
if AxC == Ax:
    print('Amplitude Constraint is Satisfied')
else:
    print('FAILED Amplitude Constraint!!')
    print('Reassigning Ax = Axc...')
    Ax = AxC

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
ax1.set(xlabel='x', ylabel='y')
ax2.plot(x, z)
ax2.set(xlabel='x', ylabel='z')
ax3.plot(y, z)
ax3.set(xlabel='y', ylabel='z')
plt.show()

# Added variables - Lindstedt-Poncairé Method

nu1 = 0
nu2 = s1*(Ax)**2+s2*(Az)**2 # REVISAR!!!!!!!!!!!
nu  = 1+nu1+nu2
tau = nu*t

tau1   = wp*tau+phi
deltam = 2-m

T      = 2*np.pi/(wp*nu)*Tconversion # Period Calculation (s)
Tdays  = T/3600/24                # s to days

# EoM - 3rd Order Richardson Expansion

xad = a21*Ax**2+a22*Az**2-Ax*np.cos(tau1)+\
    (a23*Ax**2-a24*Az**2)*np.cos(2*tau1)+(a31*Ax**3-a32*Ax*Az**2)*np.cos(3*tau1)
yad = kappa*Ax*np.sin(tau1)+\
    (b21*Ax**2-b22*Az**2)*np.sin(2*tau1)+(b31*Ax**3-b32*Ax*Az**2)*np.sin(3*tau1)
zad = deltam*Az*np.cos(tau1)+\
    deltam*d21*Ax*Az*(np.cos(2*tau1)-3)+deltam*(d32*Az*Ax**2-d31*Az**3)*np.cos(3*tau1)

x = xad*L; y = yad*L; z = zad*L

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
ax1.set(xlabel='x', ylabel='y')
ax2.plot(x, z)
ax2.set(xlabel='x', ylabel='z')
ax3.plot(y, z)
ax3.set(xlabel='y', ylabel='z')
plt.show()

# 1 orbit display

tdim = t*Tconversion
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.plot(x[tdim<T], y[tdim<T])
ax1.plot(0,0,'bx')
ax1.set(xlabel='x', ylabel='y')
ax2.plot(x[tdim<T], z[tdim<T])
ax2.plot(0,0,'bx')
ax2.set(xlabel='x', ylabel='z')
ax3.plot(y[tdim<T], z[tdim<T])
ax3.plot(0,0,'bx')
ax3.set(xlabel='y', ylabel='z')
if m == 1:
    fig.suptitle('Class I  Orbit  (L' + str(LP) +')')
else:
    fig.suptitle('Class II  Orbit  (L' + str(LP) +')')
plt.show()

if opt == 1:
    print('RCTBP: SUN-(EARTH+MOON) SYSTEM')
else:
    print('RCTBP: EARTH-MOON SYSTEM')

print('-- Results ---')
print('Ax = ' + str(AxC*L) + ' km;   Az = ' + str(AzDim) + ' km')
print('T (days) = %.5f' % Tdays)
if m == 1:
    print('Case: Northern (L' + str(LP) + ')')
else:
    print('Case: Southern (L' + str(LP) + ')')
