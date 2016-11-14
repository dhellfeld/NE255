import numpy as np
import matplotlib.pyplot as plt
import sys

dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16),])
data = np.fromfile(sys.argv[1], dtype=dt)

# access data as follows
response = np.zeros((192,3072))
for i in range(3072):
    b = np.histogram(data['DetID'][data['HPidx'] == i], 192, [0,192])
    response[:,i] += b[0]

print response
