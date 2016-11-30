# NE 255 - Homework 6 - Problem 2
# D. Hellfeld
# 12/1/16

# Imports
import numpy as np
import matplotlib.pyplot as plt
import random as random
from scipy.optimize import curve_fit

# Actual PDF(s)
def p1(x):
    return 4. * np.sqrt(1. - x**2)

def p2(x):
    return 4. * (1. / (1. + x**2))

# Simple enclosing PDF (uniform), defined on [a,b]
def g(x,a,b,const):
    if (x >= a and x <= b):
        return const
    else:
        return 0.

def G(xi,a,b,const):
    return (xi * (b-a) + a)

# Define variables
N             = np.asarray([10,50,100,500,1000,5000,10000,50000,100000])
a             = 0.
b             = 1.
const         = 4.2
pi_exact      = 3.14159
pi_estimate_1 = (np.zeros_like(N)).astype(float)
pi_estimate_2 = (np.zeros_like(N)).astype(float)
rel_error_1   = (np.zeros_like(N)).astype(float)
rel_error_2   = (np.zeros_like(N)).astype(float)

# Loop through different sample values
itr = 0
for n in N:

    # Sample
    accept_1 = 0.
    accept_2 = 0.
    for i in range(int(n)):

        xi  = random.random()
        eta = random.random()
        x   = G(xi, a, b, const)

        if (eta * g(x, a, b, const) < p1(x)):
            accept_1 += 1.
        if (eta * g(x, a, b, const) < p2(x)):
            accept_2 += 1.

    # Calculate results
    pi_estimate_1[itr] = (const*(b-a)) * (accept_1/n)
    rel_error_1[itr]   = (np.fabs(pi_estimate_1[itr] - pi_exact)/pi_exact) * 100.
    pi_estimate_2[itr] = (const*(b-a)) * (accept_2/n)
    rel_error_2[itr]   = (np.fabs(pi_estimate_2[itr] - pi_exact)/pi_exact) * 100.

    # Increment loopnum
    itr += 1

# Print results
print "Number of samples = ", N
print "Pi estimates (function 1) = ", pi_estimate_1
print "Relative errors (function 1) = ", rel_error_1, "%"
print "Pi estimates (function 2) = ", pi_estimate_2
print "Relative errors (function 2) = ", rel_error_2, "%"

# Define fit function
def fit(x,a):
    return a / np.sqrt(x)

# Fit to data
popt_1, pcov_1 = curve_fit(fit, N, rel_error_1)
popt_2, pcov_2 = curve_fit(fit, N, rel_error_2)

# Plot results with fit (function 1)
x = np.linspace(1,100000,100000)
plt.figure()
plt.plot(N, rel_error_1, linestyle='none', marker='o',label='Estimates')
plt.plot(x, fit(x,popt_1), label='$%.2f/\sqrt{N}$ fit' %popt_1)
plt.xscale('log'); plt.yscale('log')
plt.xlabel('Number of samples (N)'); plt.ylabel('Relative Error (%)')
plt.title('Relative error as function of sample number (function 1)')
plt.xlim(5,1.2e5); plt.ylim(1e-2,1e2)
plt.legend(numpoints=1)

# Plot results with fit (function 2)
plt.figure()
plt.plot(N, rel_error_2, linestyle='none', marker='o',label='Estimates')
plt.plot(x, fit(x,popt_2), label='$%.2f/\sqrt{N}$ fit' %popt_2)
plt.xscale('log'); plt.yscale('log')
plt.xlabel('Number of samples (N)'); plt.ylabel('Relative Error (%)')
plt.title('Relative error as function of sample number (function 2)')
plt.xlim(5,1.2e5); plt.ylim(1e-2,1e2)
plt.legend(numpoints=1)

# Render plots
plt.show()
