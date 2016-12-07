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

def Sequence(coinc_dets, coinc_energies):

    seqs = []
    for i in range(len(coinc_energies)):
        E1 = coinc_energies[i,0]
        E2 = coinc_energies[i,1]
        D1 = coinc_dets[i,0]
        D2 = coinc_dets[i,1]
        # Compton edge test
        if E1 <= ComptonEdge(E1+E2):
            seqs.append([E1,E2,D1,D2])
        if E2 <= ComptonEdge(E1+E2):
            seqs.append([E2,E1,D2,D1])

    return seqs


# Get the data
data = GetBinaryOutputData("../output/output_662keV_HP1200.bin")
#data = GetBinaryOutputData(sys.argv[1])
data = RemoveZeroEnergyInteractions(data)
energy  = 662
hpindex = 1200
phi,theta = np.asarray(hp.pix2ang(16,hpindex)) * (180./np.pi)

# Read in detector centers from file
detcenters = np.loadtxt('../geo/centervertices_Ring.txt')

coinc_dets     = np.array([]).reshape(0,2)
coinc_energies = np.array([]).reshape(0,2)
for i in range(len(data['EvtN']) - 1):
    if np.isclose(data['Energy'][i] + data['Energy'][i+1], energy):
        if (data['DetID'][i] != data['DetID'][i+1]):
            coinc_energies = np.vstack([coinc_energies, [data['Energy'][i], data['Energy'][i+1]]])
            coinc_dets     = np.vstack([coinc_dets,     [data['DetID'][i],  data['DetID'][i+1]]] )

sequences = Sequence(coinc_dets, coinc_energies)
#print len(sequences), "sequences"

nside = 32
[x_,y_,z_] = hp.pix2vec(nside,range(12*nside*nside))
k = zip(x_,y_,z_)

im = np.zeros(12*nside*nside)

cmap_ = plt.cm.YlGnBu_r
#cmap_ = plt.cm.jet
cmap_.set_under("w")

animate = False
if animate: plt.ion()

angunc = 3.
for i in range(len(sequences)):

    mu    = np.cos(GetScatteringAngle(sequences[i][0], sequences[i][1]))
    w     = GetConeAxis(sequences[i][2], sequences[i][3], detcenters)
    sigma = np.sin(np.arccos(mu)) * (angunc * np.pi/180.)

    val = (1. / (sigma * np.sqrt(2.*np.pi))) * np.exp(-(np.dot(k,w) - mu)**2/(2. * sigma**2))
    val[val < 1e-5] = 0
    im += val

    if animate:
        hp.cartview(im, fig=1, title="", cbar=False, cmap=cmap_)
        plt.pause(0.01)

if animate: plt.ioff()

latra = [-90,90]
lonra = [-180,180]
p = hp.cartview(im,lonra=lonra,latra=latra, return_projected_map=True)
plt.close("all")
plt.figure()
p = plt.imshow(p, cmap=cmap_, origin='lower', interpolation='nearest',extent=(lonra[0],lonra[1],latra[0],latra[1]))
#plt.scatter(theta-180, 90-phi, marker='x'); plt.xlim(lonra[0], lonra[1]); plt.ylim(latra[0], latra[1])
plt.colorbar(p, fraction=0.046, pad=0.04)
plt.title("Far Field Compton Cone Backprojection, %i keV, %i cones" %(energy, i+1))
plt.xlabel('Phi (deg)'); plt.ylabel('Theta (deg)')
plt.xticks([-180,-135,-90,-45,0,45,90,135,180]); plt.yticks([-90,-45,0,45,90])
#hp.projplot(hp.pix2ang(16,hpindex), 'k*', markersize = 8)
#hp.graticule()
plt.show()
