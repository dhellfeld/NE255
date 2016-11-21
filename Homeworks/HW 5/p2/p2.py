# Imports
import numpy as np
import matplotlib.pyplot as plt

# Within group iteration procedure (weighted diamond difference)
def InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, group):

    # Initialize new center fluxes arrays for the curret group
    psi_new = np.zeros((N,6))
    phi_new = np.zeros(N)

    # Calculate starting scalar flux (Sum 6 directions with equal weight)
    if (group == 0):
        phi = (2./np.size(mu)) * psi_group1.sum(axis=1)
    elif (group == 1):
        phi = (2./np.size(mu)) * psi_group2.sum(axis=1)
    elif (group == 2):
        phi = (2./np.size(mu)) * psi_group3.sum(axis=1)

    # Start iterating
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

            # Calculate source
            # Weight is normalized to two, and there are 18 terms (6 directions in 3 groups)
            s = q_e[group] + (2./(np.size(mu)*np.size(q_e)))  * ((sig_s[group,0]*psi_group1[i,:].sum()) + (sig_s[group,1]*psi_group2[i,:].sum()) + (sig_s[group,2]*psi_group3[i,:].sum()))

            # Calculate cetner flux
            psi_new[i,0:3] = (s + (2. * np.fabs(mu[0:3]) / (d * (1.+a))) * psi_in) / (sig_t[group] + (2. * np.fabs(mu[0:3]) / (d*(1.+a))))

            # Calcualte outgoing flux (set it equal to incoming flux for next cell)
            psi_in = ((2. / (1.+a)) * psi_new[i,0:3]) - ((1.-a)/(1.+a))*psi_in

        # Sweep in mu < 0 (left)
        for i in range(int(N)):
            s = q_e[group] + (2./(np.size(mu)*np.size(q_e)))  * ((sig_s[group,0]*psi_group1[N-i-1,:].sum()) + (sig_s[group,1]*psi_group2[N-i-1,:].sum()) + (sig_s[group,2]*psi_group3[N-i-1,:].sum()))

            psi_new[N-i-1,3:6] = (s + (2. * np.fabs(mu[3:6]) / (d * (1.+a))) * psi_in) / (sig_t[group] + (2. * np.fabs(mu[3:6]) / (d*(1.+a))))

            psi_in = ((2. / (1.+a)) * psi_new[N-i-1,3:6]) - ((1.-a)/(1.+a))*psi_in

        # Calculate scalar flux from angular flux
        # Sum 6 directions with equal weight
        phi_new = (2./np.size(mu))  * psi_new.sum(axis=1)

        # Calculate convergence criterion (l2 norm of differences)
        innercrit = np.sqrt(np.sum((phi_new - phi)**2))

        # Check convergence
        if (innercrit < 1.0e-4):
            innerconverge = True

        # Update anglular fluxes
        phi= np.copy(phi_new)
        if (group == 0):
            psi_group1 = np.copy(psi_new)
        elif (group == 1):
            psi_group2 = np.copy(psi_new)
        elif (group == 2):
            psi_group3 = np.copy(psi_new)

        # Increment inner iteration number
        inneritr += 1

    # How many inner iterations did we do?
    print '(Group = %i) Number of iterations = %i' % (group+1,inneritr)

    # Return angular and scalar flux for group
    return psi_new, phi


# -----------------------------------------------


# Define variables
a     = 0.5
mu    = np.array([0.2,0.5,0.7,-0.2,-0.5,-0.7])
sig_t = np.array([0.5,0.8,1.0])
sig_s = np.array([[0.1,0.0,0.0],[0.3,0.1,0.1],[0.1,0.3,0.3]])
q_e   = np.array([1.5,0.0,0.2])
d     = 0.1
x0    = 0.0
x1    = 2.0
N = (x1 - x0) / d
x = np.linspace((d/2.), 2.-(d/2.),N)   # center points for plotting

# Group fluxes (inital guess)
psi_group1 = np.zeros((N,6))
psi_group2 = np.zeros((N,6))
psi_group3 = np.zeros((N,6))
phi_group1 = np.zeros(N)
phi_group2 = np.zeros(N)
phi_group3 = np.zeros(N)

# Choose solver method
#method = 'GaussSeidel'
method = 'Jacobi'

# Perform outer iteration over energy groups
outerconverge = False
outeritr = 0
while (outerconverge == False):

    # Calculate current group source values
    groupsource = np.zeros(N)
    for group in range(3):
        groupsource += q_e[group] + (2./(np.size(mu)*np.size(q_e))) * (sig_s[group,0]*psi_group1.sum(axis=1)) + (sig_s[group,1]*psi_group2.sum(axis=1)) + (sig_s[group,2]*psi_group3.sum(axis=1))

    if (method == 'Jacobi'):
        # Do within group iterations (all with inital guesses)
        # Store in intermediate variable (as to not update between groups, which is Gauss-Seidel)
        psi_group1_, phi_group1_ = InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, 0)
        psi_group2_, phi_group2_ = InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, 1)
        psi_group3_, phi_group3_ = InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, 2)

        # Update fluxes
        psi_group1 = np.copy(psi_group1_); phi_group1 = np.copy(phi_group1_)
        psi_group2 = np.copy(psi_group2_); phi_group2 = np.copy(phi_group2_)
        psi_group3 = np.copy(psi_group3_); phi_group3 = np.copy(phi_group3_)

    elif (method == 'GaussSeidel'):
        # Do within group iterations, updating after each group
        psi_group1, phi_group1 = InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, 0)
        psi_group2, phi_group2 = InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, 1)
        psi_group3, phi_group3 = InnerIteration(a, mu, sig_t, sig_s, q_e, N, psi_group1, psi_group2, psi_group3, 2)

    # Calculate new groups ource values
    groupsource_new = np.zeros(N)
    for group in range(3):
        groupsource_new += q_e[group] + (2./(np.size(mu)*np.size(q_e))) * (sig_s[group,0]*psi_group1.sum(axis=1)) + (sig_s[group,1]*psi_group2.sum(axis=1)) + (sig_s[group,2]*psi_group3.sum(axis=1))

    # Calculate convergence cirterion
    outercrit = np.sqrt(np.sum((groupsource_new - groupsource)**2))

    # Check convergence
    if (outercrit < 1.0e-4):
         outerconverge = True

    # Increment outer iteration number
    outeritr += 1

    print 'Outer group iteration = %i' % outeritr


# Plot scalar fluxes
plt.figure()
plt.plot(x,phi_group1, marker='s', color='c',linestyle='none', label='Group 1')
plt.plot(x,phi_group2, marker='s', color='b',linestyle='none', label='Group 2')
plt.plot(x,phi_group3, marker='s', color='r',linestyle='none', label='Group 3')
plt.xlim(-0.1,2.1)
plt.ylim(0,8)
plt.xlabel('x'), plt.ylabel('Center $\phi$')
plt.title('Cell Center Scalar Flux Profile, $\\alpha$=%.2f' % a)
plt.legend(numpoints=1,fontsize=10)

# Render plots
plt.show()
