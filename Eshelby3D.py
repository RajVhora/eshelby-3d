import scipy.special as sp
import numpy as np
import sympy as sym

#Define constants
pi = np.pi

#Material Inputs
E = 210
nu = 0.3
mu = E/(2*(1+nu))
epsilon11 = 0.001
epsilon22 = 0.001
epsilon33 = 0.001
epsilon12 = 0
epsilon23 = 0
epsilon31 = 0

#Define 3 semi radii (condition: a1>a2>a3)
a1 = 1
a2 = 0.75
a3 = 0.5

#Find k and theta
theta = np.arcsin((1-(a3/a1)**2)**0.5)
k = ((a1**2 - a2**2)/(a1**2 - a3**2))**0.5

#Find F and E
F = sp.ellipkinc(theta,k**2)
E = sp.ellipeinc(theta,k**2)

print(theta,k)
print(F,E)

#find I1, I2, I3
I1 = 4*pi*a1*a2*a3*(F - E)/((a1**2 - a2**2)*((a1**2 - a3**2)**0.5))
I3 = 4*pi*a1*a2*a3*(((a2*(a1**2-a3**2)**0.5)/(a1*a3)) - E)/((a2**2 - a3**2)*((a1**2 - a3**2)**0.5))
I2 = 4*pi - I1 - I3
print(I1,I2,I3)

I12 = (I2-I1)/(a1**2 - a2**2)
#Finding I11 and I13
x,y = sym.symbols('x y')
eq1 = sym.Eq(3*x + I12 + y, 4*pi/a1**2)
eq2 = sym.Eq(3*(a1**2)*x + (a2**2)*I12 + (a3**2)*y, 3*I1)
sol = sym.solve([eq1, eq2], (x, y))

I11 = sol[x]
I13 = sol[y]
print(I11,I12,I13)

I23 = (I3-I2)/(a2**2 - a3**2)
#Finding I22 and I21
x,y = sym.symbols('x y')
eq1 = sym.Eq(3*x + I23 + y, 4*pi/a2**2)
eq2 = sym.Eq(3*(a2**2)*x + (a3**2)*I23 + (a1**2)*y, 3*I2)
sol = sym.solve([eq1, eq2], (x, y))

I22 = sol[x]
I21 = sol[y]
print(I21,I22,I23)

I31 = (I1-I3)/(a3**2 - a1**2)
#Finding I33 and I32
x,y = sym.symbols('x y')
eq1 = sym.Eq(3*x + I31 + y, 4*pi/a3**2)
eq2 = sym.Eq(3*(a3**2)*x + (a1**2)*I31 + (a2**2)*y, 3*I3)
sol = sym.solve([eq1, eq2], (x, y))

I33 = sol[x]
I32 = sol[y]
print(I31,I32,I33)

#Calculate sigma11
t11 = (a1**2)*(3*I11 - 3*nu*I11 + nu*I21 + nu*I31)/(8*pi*(1-nu)*(1-2*nu))
t12 = (1-2*nu)*(I1-nu*(I1+I2+I3))/(8*pi*(1-nu)*(1-2*nu))
t13 = (1-nu)/(1-2*nu)
t1 = (t11 + t12 - t13)*epsilon11

t21 = (a2**2)*(I12 - nu*I12 + 3*nu*I22 + nu*I32)/(8*pi*(1-nu)*(1-2*nu))
t22 = (1-2*nu)*(I1-nu*(I1+I2-I3))/(8*pi*(1-nu)*(1-2*nu))
t23 = (nu)/(1-2*nu)
t2 = (t21 - t22 - t23)*epsilon22

t31 = (a3**2)*(I13 - nu*I13 + nu*I23 + 3*nu*I33)/(8*pi*(1-nu)*(1-2*nu))
t32 = (1-2*nu)*(I1-nu*(I1-I2+I3))/(8*pi*(1-nu)*(1-2*nu))
t33 = (nu)/(1-2*nu)
t3 = (t31 - t32 - t33)*epsilon33

sigma11 = (t1 + t2 + t3)*2*mu

#calculate sigma22 using cyclic permutation of 1,2,3
t22 = (a2**2)*(3*I22 - 3*nu*I22 + nu*I32 + nu*I12)/(8*pi*(1-nu)*(1-2*nu))
t23 = (1-2*nu)*(I2-nu*(I2+I3+I1))/(8*pi*(1-nu)*(1-2*nu))
t21 = (1-nu)/(1-2*nu)
t2 = (t22 + t23 - t21)*epsilon22

t32 = (a3**2)*(I23 - nu*I23 + 3*nu*I33 + nu*I13)/(8*pi*(1-nu)*(1-2*nu))
t33 = (1-2*nu)*(I2-nu*(I2+I3-I1))/(8*pi*(1-nu)*(1-2*nu))
t31 = (nu)/(1-2*nu)
t3 = (t32 - t33 - t31)*epsilon33

t12 = (a1**2)*(I21 - nu*I21 + nu*I31 + 3*nu*I11)/(8*pi*(1-nu)*(1-2*nu))
t13 = (1-2*nu)*(I2-nu*(I2-I3+I1))/(8*pi*(1-nu)*(1-2*nu))
t11 = (nu)/(1-2*nu)
t1 = (t12 - t13 - t11)*epsilon11

sigma22 = (t2 + t3 + t1)*2*mu

#calculate sigma33 using cyclic permutation of 1,2,3
t33 = (a3**2)*(3*I33 - 3*nu*I33 + nu*I13 + nu*I23)/(8*pi*(1-nu)*(1-2*nu))
t31 = (1-2*nu)*(I3-nu*(I3+I1+I2))/(8*pi*(1-nu)*(1-2*nu))
t32 = (1-nu)/(1-2*nu)
t3 = (t33 + t31 - t32)*epsilon33

t13 = (a1**2)*(I31 - nu*I31 + 3*nu*I11 + nu*I21)/(8*pi*(1-nu)*(1-2*nu))
t11 = (1-2*nu)*(I3-nu*(I3+I1-I2))/(8*pi*(1-nu)*(1-2*nu))
t12 = (nu)/(1-2*nu)
t1 = (t13 - t11 - t12)*epsilon11

t23 = (a2**2)*(I32 - nu*I32 + nu*I12 + 3*nu*I22)/(8*pi*(1-nu)*(1-2*nu))
t21 = (1-2*nu)*(I3-nu*(I3-I1+I2))/(8*pi*(1-nu)*(1-2*nu))
t22 = (nu)/(1-2*nu)
t2 = (t23 - t21 - t22)*epsilon22

sigma33 = (t3 + t1 + t2)*2*mu

print(sigma11,sigma22,sigma33)

#calculate sigma12 using cyclic permutation of 1,2,3
t11 = (a1**2 + a2**2)*I12/(8*pi*(1-nu))
t12 = (1-2*nu)*(I1+I2)/(8*pi*(1-nu))
t13 = 1
sigma12 = (t11 + t12 - t13)*epsilon12*2*mu

#calculate sigma23 using cyclic permutation of 1,2,3
t22 = (a2**2 + a3**2)*I23/(8*pi*(1-nu))
t23 = (1-2*nu)*(I2+I3)/(8*pi*(1-nu))
t21 = 1
sigma23 = (t22 + t23 - t21)*epsilon23*2*mu

#calculate sigma31 using cyclic permutation of 1,2,3
t33 = (a3**2 + a1**2)*I31/(8*pi*(1-nu))
t31 = (1-2*nu)*(I3+I1)/(8*pi*(1-nu))
t32 = 1
sigma31 = (t33 + t31 - t32)*epsilon31*2*mu

print(sigma12,sigma23,sigma31)