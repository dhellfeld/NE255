import numpy as np
import matplotlib.pyplot as plt
import sys

def GetBinaryOutputData(filename):

    dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16),])
    return np.fromfile(filename, dtype=dt)

# Get the data
data = GetBinaryOutputData("../output/output_response.bin")
#data = GetBinaryOutputData(sys.argv[1])
