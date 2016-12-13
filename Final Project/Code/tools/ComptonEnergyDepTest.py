import numpy as np
import matplotlib.pyplot as plt
import sys
import healpy as hp

sequences200 = np.load('sequences200.npy')
sequences662 = np.load('sequences662.npy')


ones = np.ones(1000)

plt.figure()
plt.plot(sequences200[:1000,0]/sequences200[:1000,1], linestyle='none',marker='o', color='b', label='200 keV')
plt.plot(sequences662[:1000,0]/sequences662[:1000,1], linestyle='none',marker='o', color='r', label='662 keV')
plt.plot(ones, color='g', linewidth=2)
plt.ylim(0,3)
plt.xlabel('Sequence Number'); plt.ylabel('Ratio of First to Second Energy Deposition')
plt.title('Ratio of First Energy Deposition to Second Energy Deposition for 200 and 662 keV gamma-rays')
plt.legend(numpoints=1)
plt.show()
