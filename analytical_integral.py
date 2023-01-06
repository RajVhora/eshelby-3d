import scipy.special as sp
import numpy as np
import sympy as sym

#Define 3 semi radii
a1 = 1e20
a2 = 50.0000001
a3 = 50

#Find k and theta
theta = np.arcsin((1-(a3/a1)**2)**0.5)
k = ((a1**2 - a2**2)/(a1**2 - a3**2))**0.5

#Find F and E
F = sp.ellipkinc(theta,k)
E = sp.ellipeinc(theta,k)

#find I1, I2 and I3
