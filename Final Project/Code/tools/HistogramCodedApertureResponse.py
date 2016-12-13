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

if sys.argv[3] == "DOI":
	response = (np.histogram2d(data['DetID'] + 192.*data['DOI'], data['HPidx'], bins=(192*10,3072)))[0]
else:
	response = (np.histogram2d(data['DetID'], data['HPidx'], bins=(192,3072)))[0]

# Save response to file to avoid histogramming again
print 'Saving data'
if sys.argv[3] == "DOI":
	np.save(filename[:-4]+'_DOI', response)
else:
	np.save(filename[:-4], response)

print 'Done!'



# # Do the same for DOI
# DOI = True
# if (DOI):
#
#     # Get DOI (10 bins) system response (only full energy deposition interactions for coded aperture)
#     response_DOI10 = (np.histogram2d(data['DetID'] + 192.*data['DOI'], data['HPidx'], bins=(192*10,3072)))[0]
#
#     # Plot response
#     plt.figure()
#     im = plt.imshow(response_DOI10, origin='lower', interpolation='nearest', aspect='auto')
#     plt.colorbar(im)
#     plt.xlabel('HEALPix index (Source Location)'), plt.ylabel('Detector ID')
#     plt.title('System Response (DOI - 10 bins)')
#
#     # Pick a signal
#     #signal = response_DOI10[:,1230]
#     # Or read in input
#     signal_DOI10 = np.asarray((np.histogram(data_sig['DetID'] + 192.*data_sig['DOI'], bins=(192*10)))[0]).astype(float)
#     #signal /= signal.max()
#
#     iterations = 25
#     image_DOI10 = MLEM(response_DOI10, signal_DOI10, itr=iterations)
#     p_DOI10 = hp.cartview(image_DOI10, rot=(0,90,0), lonra=lonra,latra=latra, return_projected_map=True)
#     #plt.close("all")
#     plt.figure()
#     p_DOI10 = plt.imshow(p_DOI10, cmap=cmap_, origin='lower', interpolation='nearest',extent=(lonra[0],lonra[1],latra[0],latra[1]))
#     plt.colorbar(p, fraction=0.046, pad=0.04)
#     plt.title("Coded Aperture, %i keV, %i iterations" %(energy, iterations))
#     plt.xlabel('Phi (deg)'); plt.ylabel('Theta (deg)')
#     plt.xticks([-180,-135,-90,-45,0,45,90,135,180]); plt.yticks([-90,-45,0,45,90])
