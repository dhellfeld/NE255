import numpy as np
import matplotlib.pyplot as plt
import sys
import healpy as hp

def GetBinaryOutputData(filename):

    dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16), ('Time', np.float32)])
    return np.fromfile(filename, dtype=dt)

def RemoveZeroEnergyInteractions(data):

    return data[:][data['Energy'] != 0]

def GetScatteringAngle(E1, E2):

    return np.arccos(1. + (511./(E1+E2)) - (511./E2))

def ComptonEdge(E):
    return E / (1 + (1022/E))

def GetConeAxis(det1, det2, detcenters):

    [x,y,z] = detcenters[int(det1-1)] - detcenters[int(det2-1)]

    return [x,y,z]/np.sqrt(x**2 + y**2 + z**2)

# Get the data
data = GetBinaryOutputData("../output/output_662keV_HP912.bin")
#data = GetBinaryOutputData(sys.argv[1])
data = RemoveZeroEnergyInteractions(data)

# Read in detector centers from file
detcenters = np.loadtxt('../geo/centervertices_Ring.txt')

coinc_dets     = np.array([]).reshape(0,2)
coinc_energies = np.array([]).reshape(0,2)
for i in range(len(data['EvtN'])):
    try:
         if np.isclose(data['Energy'][i] + data['Energy'][i+1], 662):
             if (data['DetID'][i] != data['DetID'][i+1]):
                 coinc_energies = np.vstack([coinc_energies, [data['Energy'][i], data['Energy'][i+1]]])
                 coinc_dets     = np.vstack([coinc_dets,     [data['DetID'][i],  data['DetID'][i+1]]] )
    except:
        pass


def Sequence(coinc_dets, coinc_energies):

    seqs = []
    for i in range(len(coinc_energies)):
        E1 = coinc_energies[i,0]
        E2 = coinc_energies[i,1]
        D1 = coinc_dets[i,0]
        D2 = coinc_dets[i,1]
        if E1 <= ComptonEdge(E1+E2):
            seqs.append([E1,E2,D1,D2])
        if E2 <= ComptonEdge(E1+E2):
            seqs.append([E2,E1,D2,D1])

    return seqs

sequences = Sequence(coinc_dets, coinc_energies)

print len(sequences), "sequences"
nside = 16
im = np.zeros(3072)
for i in range(len(sequences)):

    angwidth = 3.

    mu    = np.cos(GetScatteringAngle(sequences[i][0], sequences[i][1]))
    w     = GetConeAxis(sequences[i][2], sequences[i][3], detcenters)
    sigma = np.sin(np.arccos(mu)) * (angwidth * np.pi/180.)

    x = hp.pix2vec(nside,range(3072))
    val = (1. / (sigma * np.sqrt(2.*np.pi))) * np.exp(-(np.dot(x,w) - mu)**2/(2. * sigma**2))

    # for pix in range(3072):
    #     x = hp.pix2vec(nside,pix)
    #
    #     val = (1. / (sigma * np.sqrt(2.*np.pi))) * np.exp(-(np.dot(x,w) - mu)**2/(2. * sigma**2))
    #     if val > 1e-5:
    #         im[pix] += val



hp.cartview(im)
plt.show()
