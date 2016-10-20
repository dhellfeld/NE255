# NE 255 - Homework 3, Problem 2c
# Daniel Hellfeld
# 10/20/16

# Imports
import numpy as np
import sys

# Function that retrieves the quadrature points, weights, and point-weight triangles
def PointsAndWeights(N):

    if N == 2:
        mu  = eta = xi = np.array([0.5773503])
        w = np.array([1.0000000])
        cells = [np.array(a) for a in [[1]]]
    elif N == 4:
        mu  = eta = xi = np.array([0.3500212, 0.8688903])
        w = np.array([0.3333333])
        cells = [np.array(a) for a in [[1], [1,1]]]
    elif N == 6:
        mu  = eta = xi = np.array([0.2666355, 0.6815076, 0.9261808])
        w = np.array([0.1761263, 0.1572071])
        cells = [np.array(a) for a in [[1], [2,2], [1,2,1]]]
    elif N == 8:
        mu  = eta = xi = np.array([0.2182179, 0.5773503, 0.7867958, 0.9511897])
        w = np.array([0.1209877, 0.0907407, 0.0925926])
        cells = [np.array(a) for a in [[1], [2,2], [2,3,2], [1,2,2,1]]]
    elif N == 12:
        mu  = eta = xi = np.array([0.1672126, 0.4595476, 0.6280191, 0.7600210, 0.8722706, 0.9716377])
        w = np.array([0.0707626, 0.0558811, 0.0373377, 0.0502819, 0.0258513])
        cells = [np.array(a) for a in [[1], [2,2], [3,4,3], [3,5,5,3], [2,4,5,4,2], [1,2,3,3,2,1]]]
    elif N == 16:
        mu = eta = xi = np.array([0.1389568, 0.3922893, 0.5370966, 0.6504264, 0.7467506, 0.8319966, 0.9092855, 0.9805009])
        w = np.array([0.0489872, 0.0413296, 0.0212326, 0.0256207, 0.0360486, 0.0144589, 0.0344958, 0.0085179])
        cells = [np.array(a) for a in [[1], [2,2], [3,5,3], [4,6,6,4], [4,7,8,7,4], [3,6,8,8,6,3], [2,5,6,7,6,5,2], [1,2,3,4,4,3,2,1]]]
    else:
        sys.exit("Error: I can only do N = 2, 4, 6, 8, 12, and 16. You entered N = %i." % N)

    return mu,eta,xi,w,cells,N

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

# Compute the quadrature integral
def QuadratureIntegral(f,mu,eta,xi,w,cells,N):

    sum_ = 0
    for octant in range(8):
        signs = GetSigns(octant)

        for i in range(N/2):
            for j in range(N/2):
                for k in range(N/2):

                    if np.isclose(mu[i]**2  + eta[j]**2 + xi[k]**2, 1):

                        val    = f(mu[i]*signs[0], eta[j]*signs[1], xi[k]*signs[2])
                        weight = w[ cells[((N/2)-1) - k][j] - 1 ]

                        sum_ += val*weight

    norm = 4.*np.pi/8.
    return norm * sum_


# --------

def Omega(mu,eta,xi):
    return np.sqrt(mu**2 + eta**2 + xi**2)

def Function1(mu,eta,xi):
    return mu + eta + xi

def Function2(mu,eta,xi):
    return mu**2


# ----------

# Choose a quadrature set (N = 2,4,6,8,12,16)
N = 8

# Calculate the integral over Omega for the three functions above
print QuadratureIntegral(Omega,*PointsAndWeights(N))       # expect to be 4pi
print QuadratureIntegral(Function1,*PointsAndWeights(N))   # expect to be 0
print QuadratureIntegral(Function2,*PointsAndWeights(N))   # expect to be 2pi*(2/3)
