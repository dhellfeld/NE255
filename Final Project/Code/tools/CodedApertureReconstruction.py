import numpy as np
import matplotlib.pyplot as plt
import sys
import healpy as hp

def MLEM(response, signal,itr = 25):

    # Remove all the empty rows from the response matrix
    #response = response[~np.all(response == 0, axis=1)]
    response[response == 0] = 1e-10
    #signal = signal[~(signal == 0.)]
    signal[signal == 0] = 1e-10

    # Get number of image bins in response
    imbins = np.shape(response)[1]

    # Remove all the empty rows from the response matrix
    response = response[~np.all(response == 0, axis=1)]

    # Normalize signal
    #signal = signal/signal.sum()

    # Initialize the image to ones
    image = np.ones(imbins)

    # Perform iterations (see Lange and Carson, 1984)
    S = np.sum(response, axis=0)  # nice to have this separate, could put directly into image line though...
    for iteration in range(1, itr + 1):
        image = (image / S) * np.dot(signal / np.dot(response, image), response)

    # Return the normalized image
    return image / np.sum(image)


# Get the system response
response = np.load(sys.argv[1])
response /= response.sum(axis=0)

# Get energy of interest
energy = float(sys.argv[2])

# Get DOI bool
if sys.argv[3] == 'DOI':
    doi = True
else:
    doi = False

# Plot response
plt.figure()
im = plt.imshow(response, origin='lower', interpolation='nearest', aspect='auto')
plt.colorbar(im)
plt.xlabel('HEALPix index (Source Location)'), plt.ylabel('Detector ID')
if doi: plt.title('System Response (DOI), %i keV' % energy)
else: plt.title('System Response (no DOI), %i keV' % energy)

if doi: response = response[:192]

# Pick a signal and get MLEM reconstruction
# Just take column of resposne
hpindex = 912
theta,phi = np.asarray(hp.pix2ang(16,hpindex-1)) * (180./np.pi)
if phi > 180: phi = -(360. - phi)
signal = response[:, hpindex]
# Or read in input
#signal = np.load(sys.argv[4])
#if doi: signal = signal[:192]
#signal /= signal.max()

iterations = 30
image = MLEM(response, signal, itr=iterations)

#cmap_ = plt.cm.YlGnBu_r
cmap_ = plt.cm.jet
cmap_.set_under("w")

latra = [-90,90]
lonra = [-180,180]
#p = hp.cartview(image, rot=(0,90,0), lonra=lonra,latra=latra, return_projected_map=True)
p = hp.cartview(image, lonra=lonra,latra=latra, return_projected_map=True)
plt.close("all")
plt.figure()
p = plt.imshow(p, cmap=cmap_, origin='lower', interpolation='nearest',extent=(lonra[1],lonra[0],latra[0],latra[1]))
plt.scatter(phi, 90-theta, marker='x'); 
plt.colorbar(p, fraction=0.046, pad=0.04)
plt.xlim(lonra[0], lonra[1]); plt.ylim(latra[0], latra[1])
plt.title("Coded Aperture, %i keV, %i iterations" %(energy, iterations))
plt.xlabel('Phi (deg)'); plt.ylabel('Theta (deg)')
plt.xticks([-180,-135,-90,-45,0,45,90,135,180]); plt.yticks([-90,-45,0,45,90])

# Render
plt.show()
