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

# Get system response (only full energy deposition interactions for coded aperture)
response = np.histogram2d(data['DetID'], data['HPidx'], bins=(192,3072))

# Plot
plt.figure()
im = plt.imshow(response[0], origin='lower', interpolation='nearest', aspect='auto')
plt.colorbar(im)
plt.show()
