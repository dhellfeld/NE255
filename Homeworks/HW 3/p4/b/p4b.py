# NE 255 - Homework 3, Problem 4b
# Daniel Hellfeld
# 10/20/2016

# Imports
import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

# Define fatcorial to return a float instead of int
def factorial_(a):
    return float(np.math.factorial(a))

# Define absolute value to return a float instead of int
def abs_(a):
    return np.fabs(a)

def Derivative(f):
    # Discrete approximation to a derivative (central difference)
    # Choose h small, but not too small to cause precision issues
    def df(x, h=1.e-4):
        return ( f(x + h) - f(x - h) ) / (2.*h)
    return df

# Function to generate the Associated Legendre Polynomials
def AssocLegendrePoly(x,l,m):
    # Terms
    a = ( ((-1.)**abs_(m)) / ( (2.**l) * factorial_(l) ) )
    b = (1. - x**2.)**(abs_(m)/2.)
    def q(y): return (y**2. - 1.)**l
    c = q

    # Take derivatives
    if (l+m > 0.):
        for i in range(int(l + abs_(m))):
            c = Derivative(c)

    # Multiply terms together
    d = a * b * c(x)

    # Check if m is negative, adjust accordingly
    if m < 0:
        p = (-1.)**abs_(m) * (factorial_(l - abs_(m)) / factorial_(l + abs_(m)))
        return p * d
    else:
        return d

def SphericalHarmonics(theta,phi,l,m):
    # Terms
    a = (-1.)**abs_(m) * np.sqrt( ((2.*l + 1.)/(4.*np.pi)) * ((factorial_(l - abs_(m)))/(factorial_(l + abs_(m)))) )
    b = AssocLegendrePoly(np.cos(theta),l,abs_(m))
    c = np.exp(1.j * abs_(m) * phi)

    # Multiply terms together
    d = a*b*c

    # Check if m is negative, adjust accordingly
    if m < 0:
        return (-1)**abs_(m) * np.conjugate(d)
    else:
        return d

# Get the signs of [mu,eta,xi] in each octant
def GetSigns(octant):

    if   octant == 0: return [+1,+1,+1]
    elif octant == 1: return [-1,+1,+1]
    elif octant == 2: return [-1,-1,+1]
    elif octant == 3: return [+1,-1,+1]
    elif octant == 4: return [+1,+1,-1]
    elif octant == 5: return [-1,+1,-1]
    elif octant == 6: return [-1,-1,-1]
    elif octant == 7: return [+1,-1,-1]


def MuEtaXiConvert(mu,eta,xi):

    theta = np.arccos(mu)
    phi = np.arctan2(xi,eta)
    return theta, phi


def D(l,m):
    if m == 0:
        deltam = 1
    else:
        deltam = 0

    return ((-1.)**m) * np.sqrt( (2.-deltam) * ((2.*l + 1.)/(4*np.pi)) * (factorial_(l-m)/factorial_(l+m)))


def Yeven(theta,phi,l,m):
    return D(l,m) * AssocLegendrePoly(np.cos(theta),l,m) * np.cos(m*phi)


def Yodd(theta,phi,l,m):
    return D(l,m) * AssocLegendrePoly(np.cos(theta),l,m) * np.sin(m*phi)


def S4_Theta_Phi():

    N = 4

    MU  = np.zeros(8*N*(N+2)/8)
    ETA = np.zeros(8*N*(N+2)/8)
    XI  = np.zeros(8*N*(N+2)/8)

    # S4 mu,eta,xi
    mu = eta = xi = np.array([0.3500212, 0.8688903])

    qq = 0
    for octant in range(8):
        signs = GetSigns(octant)
        for i in range(N/2):
            for j in range(N/2):
                for k in range(N/2):
                    if np.isclose(mu[i]**2  + eta[j]**2 + xi[k]**2, 1):
                        MU[qq]  = mu[i] * signs[0]
                        ETA[qq] = eta[i]* signs[1]
                        XI[qq]  = xi[i] * signs[2]
                        qq += 1

    return MuEtaXiConvert(MU,ETA,XI)


# ----


# Get the theta and phi for all 8 octants
[theta,phi] = S4_Theta_Phi()

for l in range(3):
    for m in range(l+1):
        print "l = ", l, ", m = ", m, " qlm = ", (4.*np.pi / 24.) * np.sum(Yeven(theta,phi,l,m))
        print "l = ", l, ", m = ", m, " slm = ", (4.*np.pi / 24.) * np.sum(Yodd(theta,phi,l,m))
