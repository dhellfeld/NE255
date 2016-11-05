# python python_stitcher.py <nside> <doibins> <outputfiles>

import numpy as np
import matplotlib.pyplot as plt
import pickle as pickle
import sys as sys

nside = int(sys.argv[1])
hpi = 12*nside*nside

doibins = int(sys.argv[2])

response = np.zeros(192*doibins*hpi)

for i in range(len(sys.argv)-3):
    
    fName = sys.argv[i+3]
    response += np.fromfile(fName, dtype=np.dtype(np.uint16))

response = response.reshape(192*doibins,hpi)


# dump response to file
#pickle.dump(response, open("systemresponse.p","wb"))     (slow for large array!)
np.save("output/systemresponse",response)



