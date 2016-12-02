import numpy as np
import matplotlib.pyplot as plt
import sys
import healpy as hp

def GetBinaryOutputData(filename):

    dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16), ('Time', np.float32)])
    return np.fromfile(filename, dtype=dt)

def FullEnergyAbsorptions(data, fullenergy):

    return data[:][data['Energy'] == fullenergy]

def MLEM(response, signal,itr = 25):

    # Get number of image bins in response
    imbins = np.shape(response)[1]

    # Remove all the empty rows from the response matrix
    response = response[~np.all(response == 0, axis=1)]

    # Normalize signal
    signal = signal/signal.sum()

    # Initialize the image to ones
    image = np.ones(imbins)

    # Perform iterations (see Lange and Carson, 1984)
    S = np.sum(response, axis=0)  # nice to have this separate, could put directly into image line though...
    for iteration in range(1, itr + 1):
        image = (image / S) * np.dot(signal / np.dot(response, image), response)

    # Return the normalized image
    return image / np.sum(image)


# Get the data
data = GetBinaryOutputData("../output/output_response.bin")
#data = GetBinaryOutputData(sys.argv[1])

# Pull out only full energy absorptions (60 keV)
data = FullEnergyAbsorptions(data, 60.0)

# Get non-DOI system response (only full energy deposition interactions for coded aperture)
response_noDOI = (np.histogram2d(data['DetID'], data['HPidx'], bins=(192,3072)))[0]

# Plot response
plt.figure()
im = plt.imshow(response_noDOI, origin='lower', interpolation='nearest', aspect='auto')
plt.colorbar(im)
plt.xlabel('HEALPix index (Source Location)'), plt.ylabel('Detector ID')
plt.title('System Response (no DOI)')

# Pick a signal and get MLEM reconstruction
signal = response_noDOI[:,1230]
image = MLEM(response_noDOI, signal)
hp.mollview(image)

# Do the same for DOI
DOI = False
if (DOI):

    # Get DOI (10 bins) system response (only full energy deposition interactions for coded aperture)
    response_DOI10 = (np.histogram2d(data['DetID'] + 192.*data['DOI'], data['HPidx'], bins=(192*10,3072)))[0]

    # Plot response
    plt.figure()
    im = plt.imshow(response_DOI10, origin='lower', interpolation='nearest', aspect='auto')
    plt.colorbar(im)
    plt.xlabel('HEALPix index (Source Location)'), plt.ylabel('Detector ID')
    plt.title('System Response (DOI - 10 bins)')

    # Pick a signal and get MLEM reconstruction
    signal = response_DOI10[:,1230]
    image = MLEM(response_DOI10, signal)
    hp.mollview(image)

# Render
plt.show()
