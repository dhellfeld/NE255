# NE 255 - Homework 4, Problem 3c-d
# Daniel Hellfeld
# 11/8/16

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Define variables
a = 0.0 # -0.5, 0.0, 0.5
mu = np.asarray([0.2,0.7,-0.2,-0.7])
sig_t = 1.0
sig_s = 0.9
q = 1.0
d_array = np.asarray([0.4,0.2,0.125,0.1,0.08])
weight = np.asarray([0.5,0.5,0.5,0.5])   # (made up quaderature set, equiprobable, normalized to 2)

# Loop through mesh spacings
loopnum = 1   # for subplots
for d in d_array:

    # Define more variables
    N = 2. / d
    x = np.linspace((d/2.), 2.-(d/2.),N)   # get center points for plotting

    # Initialize center fluxes arrays (2D now b/c two angles)
    psi_cent_right = np.zeros((N,2))
    psi_cent_left  = np.zeros((N,2))
    psi_cent_right_new = np.zeros((N,2))
    psi_cent_left_new  = np.zeros((N,2))
    phi = np.zeros(N)
    phi_new = np.zeros(N)

    # Iterate
    converge = False
    itr = 0
    while (converge == False):

        # Boundary condition
        psi_in = np.asarray([2.,2.])
        psi_out = np.asarray([0.,0.])

        # Sweep in mu > 0 (right)
        for i in range(int(N)):
            s = q + sig_s * (weight[0]*psi_cent_right[i,0] + weight[1]*psi_cent_right[i,1] + weight[2]*psi_cent_left[i,0] +  weight[3]*psi_cent_left[i,1])

            psi_cent_right_new[i,0] = (s + (2. * np.fabs(mu[0]) / (d * (1.+a))) * psi_in[0]) / (sig_t + (2. * np.fabs(mu[0]) / (d*(1.+a))))
            psi_cent_right_new[i,1] = (s + (2. * np.fabs(mu[1]) / (d * (1.+a))) * psi_in[1]) / (sig_t + (2. * np.fabs(mu[1]) / (d*(1.+a))))

            psi_out[0] = ((2. / (1.+a)) * psi_cent_right[i,0]) - ((1.-a)/(1.+a))*psi_in[0]
            psi_out[1] = ((2. / (1.+a)) * psi_cent_right[i,1]) - ((1.-a)/(1.+a))*psi_in[1]

            psi_in[0] = psi_out[0]
            psi_in[1] = psi_out[1]

        # Sweep in mu < 0 (left)
        for i in range(int(N)):
            s = q + sig_s * (weight[0]*psi_cent_right[N-i-1,0] + weight[1]*psi_cent_right[N-i-1,1] + weight[2]*psi_cent_left[N-i-1,0] +  weight[3]*psi_cent_left[N-i-1,1])

            psi_cent_left_new[N-i-1,0] = (s + (2. * np.fabs(mu[2]) / (d*(1.-a))) * psi_in[0]) / (sig_t + (2. * np.fabs(mu[2]) / (d*(1.-a))))
            psi_cent_left_new[N-i-1,1] = (s + (2. * np.fabs(mu[3]) / (d*(1.-a))) * psi_in[1]) / (sig_t + (2. * np.fabs(mu[3]) / (d*(1.-a))))

            psi_out[0] = ((2./(1.-a)) * psi_cent_left[N-i-1,0]) - ((1.+a)/(1.-a))*psi_in[0]
            psi_out[1] = ((2./(1.-a)) * psi_cent_left[N-i-1,1]) - ((1.+a)/(1.-a))*psi_in[1]

            psi_in[0] = psi_out[0]
            psi_in[1] = psi_out[1]

        # Calculate scalar flux from angular flux (quaderature)
        phi_new = weight[0]*psi_cent_right_new[:,0] + weight[1]*psi_cent_right_new[:,1] + weight[2]*psi_cent_left_new[:,0] + weight[3]*psi_cent_left_new[:,1]

        # Calculate convergence criterion
        crit = np.sqrt(np.sum((phi_new - phi)**2))

        # Update fluxes
        psi_cent_right = psi_cent_right_new
        psi_cent_left = psi_cent_left_new
        phi = phi_new

        # Check convergence
        if (crit < 0.001):
            converge = True
        else:
            itr += 1

    # How many iterations did we do?
    print '(mesh spacing = %.3f) Number of iterations = %i' % (d,itr)

    # Plot angular flux
    plt.subplot(5,2,2*loopnum-1)
    plt.subplots_adjust(hspace = 0.0, wspace=0.2)
    plt.plot(x,psi_cent_right[:,0], marker="o", linestyle='none',label='$\mu=0.2$')
    plt.plot(x,psi_cent_right[:,1], marker="o", linestyle='none', label='$\mu=0.7$')
    plt.plot(x,psi_cent_left[:,0], marker="o", linestyle='none',label='$\mu=-0.2$')
    plt.plot(x,psi_cent_left[:,1], marker="o", linestyle='none', label='$\mu=-0.7$')
    plt.xlim(-0.1,2.4)
    plt.ylim(np.min(psi_cent_left)-0.5*np.max(psi_cent_left),np.max(psi_cent_right)+0.3*np.max(psi_cent_right))
    plt.xlabel('x'), plt.ylabel('Center $\psi$, ($\Delta_i$=%.2f)\n(Iterations = %i)' % (d,itr))
    if (loopnum == 1): plt.title('Cell Center Angular Flux Profile, $\\alpha$=%0.2f' % a)
    plt.legend(numpoints=1,fontsize=10)

    # Plot scalar flux
    plt.subplot(5,2,2*loopnum)
    plt.plot(x,phi, marker='s', color='c',linestyle='none')
    plt.xlim(-0.1,2.1)
    plt.ylim(np.min(phi)-0.3*np.max(phi),np.max(phi)+0.3*np.max(phi))
    plt.xlabel('x'), plt.ylabel('Center $\phi$, ($\Delta_i$=%.2f)\n(Iterations = %i)' % (d,itr))
    if (loopnum == 1): plt.title('Cell Center Scalar Flux Profile, $\\alpha$=%.2f' % a)

    # Increment loopnum
    loopnum += 1

# Render plots
plt.show()
