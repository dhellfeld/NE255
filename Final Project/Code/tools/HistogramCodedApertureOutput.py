import numpy as np
import matplotlib.pyplot as plt
import sys

def GetBinaryOutputData(filename):

    dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16), ('Time', np.float32)])
    return np.fromfile(filename, dtype=dt)

def FullEnergyAbsorptions(data, fullenergy):

    return data[:][data['Energy'] == fullenergy]

# Get the data
filename = sys.argv[1]
print 'Parsing output file'
data = GetBinaryOutputData(filename)

# Pull out only full energy absorptions
energy = float(sys.argv[2])
print 'Pulling out full energy absorptions'
data = FullEnergyAbsorptions(data, energy)

# Get non-DOI system response (only full energy deposition interactions for coded aperture)
print 'Histogramming data'
output = np.asarray((np.histogram(data_sig['DetID'], bins=(192)))[0]).astype(float)

# Save response to file to avoid histogramming again
print 'Saving data'
np.save(filename[:-4], output)

print 'Done!'
