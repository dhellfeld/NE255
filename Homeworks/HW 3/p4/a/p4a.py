# NE 255 - Homework 3, Problem 4a
# Daniel Hellfeld
# 10/20/16

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
    def df(x, h=1.e-2):
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

def SH_plot(l, m, space='real'):

    # Define our theta and phi spacing
    theta = np.linspace(0, np.pi, 100.)
    phi = np.linspace(0, 2.*np.pi, 100.)

    # Make into a mesh
    [THETA,PHI] = np.meshgrid(theta,phi)

    # Define x,y,z from our theta and phi
    r = 1.
    x = r*np.sin(THETA)*np.cos(PHI)
    y = r*np.sin(THETA)*np.sin(PHI)
    z = r*np.cos(THETA)

    # Determine whether we want real or imaginary part (real is default)
    if space == "real":
        ww = np.real(SphericalHarmonics(THETA,PHI,l,m))
    elif space == 'imag':
        if m == 0: return  # imaginary part is zero everywhere when m = 0
        else: ww = np.imag(SphericalHarmonics(THETA,PHI,l,m))
    else:
        print "I do not recognize your input space command (either 'real' or 'imag')"


    # Get the absolute value of the SH functions
    # This is used to generate the distored sphere images
    qq = abs(ww)

    # Plot distored sphere with intensity
    mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))
    p = mlab.mesh(x*qq, y*qq, z*qq, scalars=ww, colormap='jet')
    mlab.colorbar(p, orientation="vertical", title="%s(Y$%i%i$)" %(space,l,m))
    mlab.axes(extent=[-0.6, 0.6, -0.6, 0.6, -0.6, 0.6], nb_labels=3)

    # Plot a scaled down unit sphere with intensity
    scale = 0.3
    p1 = mlab.mesh(x*scale, y*scale, z*scale + 1.5, scalars=ww, colormap='jet')

    mlab.view(45, 70, 7)



# ------


l_array = np.array((0,1,2))

# Plot the real and imaginary values
for l in l_array:
    if l == 0 :
        m = 0xc
        SH_plot(l,m,'real')
        SH_plot(l,m,'imag')
    else:
        for m in range(-l, l+1):
            SH_plot(l,m,'real')
            SH_plot(l,m,'imag')

mlab.show()
