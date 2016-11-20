# NE 255 - Homework 4, Problem 2
# Daniel Hellfeld
# 11/22/16

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Define variables
a     = 0.5
mu    = np.array([0.2,0.5,0.7,-0.2,-0.5,-0.7])
sig_t = np.array([0.5,0.8,1.0])
sig_s = np.array([[0.1,0.0,0.0],[0.3,0.1,0.1],[0.1,0.3,0.3]])
q_e   = np.array([1.5,0.0,0.2])
d     = 0.1
x0    = 0.0
x1    = 2.0

# Define more variables
N = (x1 - x0) / d
x = np.linspace((d/2.), 2.-(d/2.),N)   # get center points for plotting

# Group fluxes
psi_group1 = np.zeros((N,6))
psi_group2 = np.zeros((N,6))
psi_group3 = np.zeros((N,6))
phi_group1 = np.zeros(N)
phi_group2 = np.zeros(N)
phi_group3 = np.zeros(N)

outerconverge = False
outeritr = 0
while (outerconverge == False):

    groupsource = np.zeros(N)
    for group in range(3):
        groupsource += q_e[group] + (2./6.) * (sig_s[group,0]*psi_group1.sum(axis=1)) + (sig_s[group,1]*psi_group2.sum(axis=1)) + (sig_s[group,2]*psi_group3.sum(axis=1))

    for group in range(3):

        # Initialize center fluxes arrays
        psi_new = np.zeros((N,6))
        phi_new = np.zeros(N)

        # Iterate
        innerconverge = False
        inneritr = 0
        while (innerconverge == False):

            # Boundary condition
            if (group == 0):
                psi_in  = np.array([0.5,0.5,0.5])
            else:
                psi_in  = np.array([0.0,0.0,0.0])

            # Sweep in mu > 0 (right)
            for i in range(int(N)):
                s = q_e[group] + (2./6.) * ((sig_s[group,0]*psi_group1[i,:].sum()) + (sig_s[group,1]*psi_group2[i,:].sum()) + (sig_s[group,2]*psi_group3[i,:].sum()))

                psi_new[i,0:3] = (s + (2. * np.fabs(mu[0:3]) / (d * (1.+a))) * psi_in) / (sig_t[group] + (2. * np.fabs(mu[0:3]) / (d*(1.+a))))

                psi_in = ((2. / (1.+a)) * psi_new[i,0:3]) - ((1.-a)/(1.+a))*psi_in


            # Sweep in mu < 0 (left)
            for i in range(int(N)):
                s = q_e[group] + (2./6.) * ((sig_s[group,0]*psi_group1[N-i-1,:].sum()) + (sig_s[group,1]*psi_group2[N-i-1,:].sum()) + (sig_s[group,2]*psi_group3[N-i-1,:].sum()))

                psi_new[N-i-1,3:6] = (s + (2. * np.fabs(mu[3:6]) / (d * (1.+a))) * psi_in) / (sig_t[group] + (2. * np.fabs(mu[3:6]) / (d*(1.+a))))

                psi_in = ((2. / (1.+a)) * psi_new[N-i-1,3:6]) - ((1.-a)/(1.+a))*psi_in


            # Calculate scalar flux from angular flux (quaderature)
            phi_new = (2./6.) * psi_new.sum(axis=1)

            # Calculate convergence criterion (l2 norm of differences)
            if (group == 0):
                innercrit = np.sqrt(np.sum((phi_new - phi_group1)**2))
                psi_group1 = np.copy(psi_new)
                phi_group1 = np.copy(phi_new)
            elif (group == 1):
                innercrit = np.sqrt(np.sum((phi_new - phi_group2)**2))
                psi_group2 = np.copy(psi_new)
                phi_group2 = np.copy(phi_new)
            elif (group == 2):
                innercrit = np.sqrt(np.sum((phi_new - phi_group3)**2))
                psi_group3 = np.copy(psi_new)
                phi_group3 = np.copy(phi_new)

            # Check convergence
            if (innercrit < 0.0001):
                innerconverge = True

            # Increment iteration number
            inneritr += 1

        # How many iterations did we do?
        print '(Group = %i) Number of iterations = %i' % (group+1,inneritr)

    groupsource_new = np.zeros(N)
    for group in range(3):
        groupsource_new += q_e[group] + (2./6.) * (sig_s[group,0]*psi_group1.sum(axis=1)) + (sig_s[group,1]*psi_group2.sum(axis=1)) + (sig_s[group,2]*psi_group3.sum(axis=1))

    outercrit = np.sqrt(np.sum((groupsource_new - groupsource)**2))

    if (outercrit < 1.0e-4):
         outerconverge = True

    outeritr += 1

    print 'Outer group interation = %i' % outeritr

# Plot scalar flux
plt.figure()
plt.plot(x,phi_group1, marker='s', color='c',linestyle='none', label='Group 1')
plt.plot(x,phi_group2, marker='s', color='b',linestyle='none', label='Group 2')
plt.plot(x,phi_group3, marker='s', color='r',linestyle='none', label='Group 3')
plt.xlim(-0.1,2.1)
plt.ylim(0,25)
plt.xlabel('x'), plt.ylabel('Center $\phi$')
plt.title('Cell Center Scalar Flux Profile, $\\alpha$=%.2f' % a)
plt.legend(numpoints=1,fontsize=10)

# Render plots
plt.show()
