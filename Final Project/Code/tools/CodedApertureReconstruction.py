import numpy as np
import matplotlib.pyplot as plt
import sys

def GetBinaryOutputData(filename):

    dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16),])
    return np.fromfile(filename, dtype=dt)

def FullEnergyAbsorptions(data, fullenergy):

    return data[:][data['Energy'] == fullenergy]

# Get the data
data = GetBinaryOutputData("../output/output_response.bin")
#data = GetBinaryOutputData(sys.argv[1])

# Pull out only full energy absorptions (60 keV)
data = FullEnergyAbsorptions(data, 60.0)

# Get non-DOI system response (only full energy deposition interactions for coded aperture)
response_noDOI = (np.histogram2d(data['DetID'], data['HPidx'], bins=(192,3072)))[0]

# Plot
plt.figure()
im = plt.imshow(response_noDOI, origin='lower', interpolation='nearest', aspect='auto')
plt.colorbar(im)
plt.xlabel('HEALPix index (Source Location)'), plt.ylabel('Detector ID')
plt.title('System Response (no DOI)')

# Do the same for DOI
DOI = True
if (DOI):

    # Get DOI (10 bins) system response (only full energy deposition interactions for coded aperture)
    response_DOI10 = (np.histogram2d(data['DetID'] + 192.*data['DOI'], data['HPidx'], bins=(192*10,3072)))[0]

    # Plot
    plt.figure()
    im = plt.imshow(response_DOI10, origin='lower', interpolation='nearest', aspect='auto')
    plt.colorbar(im)
    plt.xlabel('HEALPix index (Source Location)'), plt.ylabel('Detector ID')
    plt.title('System Response (DOI - 10 bins)')

# Render
plt.show()
