import numpy as np
import matplotlib.pyplot as plt

def KleinNishina(E,theta):
    a = 1./137.04
    r = 0.38616e-13
    P = 1. / (1. + (E/511.)*(1-np.cos(theta)))

    return ((a**2 * r**2 * P**2) / 2.) * (P + (1./P) - 1. + np.cos(theta)**2)


theta = np.linspace(0, 2.*np.pi, 10000)
ax = plt.subplot(111, projection='polar')
ax.plot(theta, KleinNishina(1,    theta), color='black',  linewidth=3, label='1 keV')
ax.plot(theta, KleinNishina(60,   theta), color='orange', linewidth=3, label='60 keV')
ax.plot(theta, KleinNishina(200,  theta), color='g',      linewidth=3, label='200 keV')
ax.plot(theta, KleinNishina(662,  theta), color='b',      linewidth=3, label='662 keV')
ax.plot(theta, KleinNishina(1460, theta), color='r',      linewidth=3, label='1.46 MeV')
ax.plot(theta, KleinNishina(3000, theta), color='y',      linewidth=3, label='3 MeV')
#ax.set_rmax(8e-32)
ax.grid(True)
plt.legend()
plt.show()
