# python plotter.py <filename.npy>

import numpy as np
import matplotlib.pyplot as plt
import sys as sys

fName = sys.argv[1]

# load in response
response = np.load(fName)

dets = np.shape(response)[0]
hpi  = np.shape(response)[1]

# Plot the response
plt.figure()
im = plt.imshow(response,interpolation="nearest", origin="lower", aspect="auto")
plt.colorbar(im)
plt.xlabel("HP index"); plt.ylabel("Detector ID"); plt.title("System Reponse - 10 DOI bins")
plt.show()


if dets == 1920:

	# Collapse to 5 depths
	response_DOI5 = np.zeros((192*5,hpi))
	
	for i in range(192):
	    response_DOI5[(192*0)+i,:] = response[i+192*0,:] + response[i+192*1,:]
	    response_DOI5[(192*1)+i,:] = response[i+192*2,:] + response[i+192*3,:]
	    response_DOI5[(192*2)+i,:] = response[i+192*4,:] + response[i+192*5,:]
	    response_DOI5[(192*3)+i,:] = response[i+192*6,:] + response[i+192*7,:]
	    response_DOI5[(192*4)+i,:] = response[i+192*8,:] + response[i+192*9,:]
	
	# Plot the response
	plt.figure()
	im = plt.imshow(response_DOI5,interpolation="nearest", origin="lower", aspect="auto")
	plt.colorbar(im)
	plt.xlabel("HP index"); plt.ylabel("Detector ID"); plt.title("System Reponse - 5 DOI bins")
	plt.show()
	
	
	
	# Collapse to 2 depths (I+O)
	response_DOI2 = np.zeros((192*2,hpi))
	
	for i in range(192):
	    response_DOI2[i,:]     = response[i+192*0,:] + response[i+192*1,:] + response[i+192*2,:] + response[i+192*3,:] + response[i+192*4,:]
	    response_DOI2[192+i,:] = response[i+192*5,:] + response[i+192*6,:] + response[i+192*7,:] + response[i+192*8,:] + response[i+192*9,:]
	
	# Plot the response
	plt.figure()
	im = plt.imshow(response_DOI2,interpolation="nearest", origin="lower", aspect="auto")
	plt.colorbar(im)
	plt.xlabel("HP index"); plt.ylabel("Detector ID"); plt.title("System Reponse - 2 DOI bins")
	plt.show()
	
	
	
	# Collapse to no DOI
	response_noDOI = np.zeros((192,hpi))
	
	for i in range(192):
	    response_noDOI[i,:] = response[i,:] + response[i+192*1,:] + response[i+192*2,:] + response[i+192*3,:] + response[i+192*4,:] + response[i+192*5,:] + response[i+192*6,:] + response[i+192*7,:] + response[i+192*8,:] + response[i+192*9,:]
	
	# Plot the response
	plt.figure()
	im = plt.imshow(response_noDOI,interpolation="nearest", origin="lower", aspect="auto")
	plt.colorbar(im)
	plt.xlabel("HP index"); plt.ylabel("Detector ID"); plt.title("System Reponse - no DOI")
	plt.show()
	