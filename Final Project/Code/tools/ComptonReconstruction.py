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
    return E - (E / (1. + (2.*E/511.)))

def GetConeAxis(det1, det2, detcenters):

    [x,y,z] = detcenters[int(det1-1)] - detcenters[int(det2-1)]

    return [x,y,z]/np.sqrt(x**2 + y**2 + z**2)

def KleinNishina(E,theta):
    #a = 1./137.04
    #r = 0.38616e-13
    a = 1.
    r = 1.
    P = 1. / (1. + (E/511.)*(1.-np.cos(theta)))

    return ((a**2 * r**2 * P**2) / 2.) * (P + (1./P) - 1. + np.cos(theta)**2)

def GetSequences(data):

    sequences = []
    for i in range(len(data['EvtN']) - 1):
        E1 = data['Energy'][i]
        E2 = data['Energy'][i+1]
        D1 = data['DetID'][i]
        D2 = data['DetID'][i+1]
        E_ce = ComptonEdge(E1+E2)
        if np.isclose(E1 + E2, energy):
            if (D1 != D2):
                if (E1 <= E_ce) and (E2 <= E_ce):
                    # 1 means ambiguous track
                    sequences.append([E1,E2,D1,D2,1])
                    sequences.append([E2,E1,D2,D1,1])
                elif E1 <= E_ce:
                    # zero means unambiguous track
                    sequences.append([E1,E2,D1,D2,0])
                elif E2 <= E_ce:
                    # zero means unambiguous track
                    sequences.append([E2,E1,D2,D1,0])

    return sequences

def LeverArm(D1, D2, detcenters):

    [x,y,z] = detcenters[int(D1-1)] - detcenters[int(D2-1)]

    return np.sqrt(x**2 + y**2 + z**2)


# Get the data
data = GetBinaryOutputData("../output/output_662kev_full_HP1575.bin")
#data = GetBinaryOutputData(sys.argv[1])
data = RemoveZeroEnergyInteractions(data)
energy  = 662
hpindex = 1575
theta,phi = np.asarray(hp.pix2ang(16,hpindex-1)) * (180./np.pi)
if phi > 180: phi = -(360. - phi)

# Read in detector centers from file
detcenters = np.loadtxt('../geo/centervertices_Ring.txt')

# Get sequences
sequences = GetSequences(data)

nside = 32
[x_,y_,z_] = hp.pix2vec(nside,range(12*nside*nside))
k = zip(x_,y_,z_)

im = np.zeros(12*nside*nside)

#cmap_ = plt.cm.YlGnBu_r
cmap_ = plt.cm.jet
cmap_.set_under("w")

animate = True
if animate: plt.ion()

angunc = 3.
for i in range(50):
    i += 1000

    mu     = np.cos(GetScatteringAngle(sequences[i][0], sequences[i][1]))
    w      = GetConeAxis(sequences[i][2], sequences[i][3], detcenters)
    sigma  = np.sin(np.arccos(mu)) * (angunc * np.pi/180.)
    L      = LeverArm(sequences[i][2], sequences[i][3], detcenters) / 10.
    weight = 1.

    # If track is ambiguous, apply a weighting proportional to KN probability
    if sequences[i][4] == 1:
        a = KleinNishina(sequences[i][0]+sequences[i][1], mu)
        b = KleinNishina(sequences[i][0]+sequences[i][1], np.cos(GetScatteringAngle(sequences[i][1], sequences[i][0])))
        weight *= a/(a+b)

    # Lever arm waiting (longer lever arm = higher waiting)
    weight *= L**2

    val = (weight / (sigma * np.sqrt(2.*np.pi))) * np.exp(-(np.dot(k,w) - mu)**2/(2. * sigma**2))
    val[val < 1e-5] = 0
    im += val

    if animate:
        hp.cartview(im, fig=1, title="", cbar=False, cmap=cmap_)
        plt.pause(0.01)

if animate: plt.ioff()

latra = [-90,90]
lonra = [-180,180]
p = hp.cartview(im,lonra=lonra,latra=latra, return_projected_map=True)
hp.projplot(hp.pix2ang(16,hpindex-1), 'k*', markersize = 8)
#hp.graticule()
plt.close("all")
plt.figure()
p = plt.imshow(p, cmap=cmap_, origin='lower', interpolation='nearest', extent=(lonra[1],lonra[0],latra[0],latra[1]))
plt.scatter(phi, 90-theta, marker='x'); plt.xlim(lonra[0], lonra[1]); plt.ylim(latra[0], latra[1])
plt.colorbar(p, fraction=0.046, pad=0.04)
plt.title("Far Field Compton Cone Backprojection, %i keV, %i cones" %(energy, i+1))
plt.xlabel('Phi (deg)'); plt.ylabel('Theta (deg)')
plt.xticks([-180,-135,-90,-45,0,45,90,135,180]); plt.yticks([-90,-45,0,45,90])
plt.show()
