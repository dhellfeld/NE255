import numpy as np
import matplotlib.pyplot as plt
import sys

def GetBinaryOutputData(filename):

    dt = np.dtype([('EvtN', np.uint32), ('HitN', np.uint8), ('TrackID', np.uint8), ('Energy', np.float32), ('DetID', np.uint8), ('Proc', np.uint8), ('DOI', np.uint8), ('HPidx', np.uint16), ('Time', np.float32)])
    return np.fromfile(filename, dtype=dt)

def RemoveZeroEnergyInteractions(data):

    return data[:][data['Energy'] != 0]

def GetScatteringAngle(E1, E2):

    #return np.arccos(1. - (511.*E1)/((E1+E2) * E2)) * 180./np.pi
    return np.arccos(1. + (511./(E1+E2)) - (511./E2))* 180./np.pi


# Get the data
data = GetBinaryOutputData("../output/output_662keV_HP45.bin")
#data = GetBinaryOutputData(sys.argv[1])
data = RemoveZeroEnergyInteractions(data)

# Read in detector centers from file
detcenters = np.loadtxt('../geo/centervertices_Ring.txt')
# Detector number 3, y position detcenters[2,1]

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

print coinc_energies[0,0], coinc_energies[0,1]
print GetScatteringAngle(coinc_energies[0,0], coinc_energies[0,1])
print GetScatteringAngle(coinc_energies[0,1], coinc_energies[0,0])
