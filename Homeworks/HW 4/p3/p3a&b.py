# NE 255 - Homework 4, Problem 3a-b
# Daniel Hellfeld
# 11/8/16

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Define variables
a = 0.0     # -0.9, -0.5, 0, 0.25, 0.5, 0.9
mu = np.asarray([0.1,-0.1])
sig_t = 1.0
d_array = np.asarray([0.4,0.2,0.125,0.1,0.08])
weight = np.asarray([1.,1.])   # (made up quaderature set, equiprobable, normalized to 2)

# Loop through mesh spacings
loopnum = 1   # for subplots
for d in d_array:

    # Define more variables
    N = 2. / d
    x = np.linspace((d/2.), 2.-(d/2.),N)   # get center points for plotting

    # Initialize center fluxes arrays
    psi_cent_right = np.zeros(N)
    psi_cent_left  = np.zeros(N)

    # Boundary condition
    psi_in = 2.

    # Sweep in mu > 0 (right)
    for i in range(int(N)):
        psi_cent_right[i] = ((2. * np.abs(mu[0]) / (d * (1.+a))) * psi_in) / (sig_t + (2. * np.abs(mu[0]) / (d*(1.+a))))
        psi_out = ((2. / (1.+a)) * psi_cent_right[i]) - ((1.-a)/(1.+a))*psi_in
        psi_in = psi_out

    # Sweep in mu < 0 (left)
    for i in range(int(N)):
        psi_cent_left[N-i-1] = ((2. * np.abs(mu[1]) / (d*(1.-a))) * psi_in) / (sig_t + (2. * np.abs(mu[1]) / (d*(1.-a))))
        psi_out = ((2./(1.-a)) * psi_cent_left[N-i-1]) - ((1.+a)/(1.-a))*psi_in
        psi_in = psi_out

    # Calculate scalar flux from angular flux (quaderature)
    phi = weight[0]*psi_cent_left + weight[1]*psi_cent_right

    # Plot angular flux
    plt.subplot(5,2,2*loopnum-1)
    plt.subplots_adjust(hspace = 0.0, wspace=0.2)
    plt.plot(x,psi_cent_right, marker="o", linestyle='none',label='$\mu>0$')
    plt.plot(x,psi_cent_left, marker="o",color='r', linestyle='none', label='$\mu<0$')
    plt.xlim(-0.1,2.1)
    plt.ylim(np.min(phi)-0.5*np.max(phi),np.max(psi_cent_right)+0.3*np.max(psi_cent_right))
    plt.xlabel('x'), plt.ylabel('Center $\psi$, ($\Delta_i$=%.2f)' % d)
    if (loopnum == 1): plt.title('Cell Center Angular Flux Profile, $\\alpha$=%0.2f' % a)
    plt.legend(numpoints=1,fontsize=10)

    # Plot scalar flux
    plt.subplot(5,2,2*loopnum)
    plt.plot(x,phi, marker='s', color='c',linestyle='none')
    plt.xlim(-0.1,2.1)
    plt.ylim(np.min(phi)-0.3*np.max(phi),np.max(phi)+0.3*np.max(phi))
    plt.xlabel('x'), plt.ylabel('Center $\phi$, ($\Delta_i$=%.2f)' % d)
    if (loopnum == 1): plt.title('Cell Center Scalar Flux Profile, $\\alpha$=%.2f' % a)

    # Increment loopnum
    loopnum += 1

# Render plots
plt.show()
