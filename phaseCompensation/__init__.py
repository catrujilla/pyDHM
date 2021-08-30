# libraries
import numpy as np
from math import pi

def info():
    
    print ("The available phase compensation algorithms are Full ROI Search (FRS), Efficient ROI Search (ERS) and the Cost Function ROI (CFR)")
    
    return
    
# Function to calcule the intensity representation of a given complex field
# Inputs:
# inp - The input complex field
# log - boolean variable to determine if a log representation is applied
def intensity (inp, log):
    
    out = np.abs(inp)
    out = out*out
    
    if log == True:
        out = 20 * np.log(out)
	
    return out

#Spatial filtering process
def spatialFiltering (inp, M, N, upper):
    
    field_spec = np.fft.fftshift(inp)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)
    
    pow_espec = intensity(field_spec, True)  
    
    #Defining width and hieght of point of maxima in the ROI determination
    #10% is estimated to be the proper distance to separated an in-line hologram from an off-axis one in the Spectrum
    M_search = int(M-(M*0.1))
    N_search = int(N-(N*0.1))
    #print (M_search)
    
    #Location and size of the ROI for the +1 or -1 diffraction order info
    
    #Location of the center of the ROI in the DHM hologram spectrum
    if (upper == True):
        #Find the location of the maxima point in the upper part of the DHM hologram spectrum
        indI = np.unravel_index(np.argmax(pow_espec[1:(int(M/2)-M_search), 1:N], axis=None), pow_espec[1:(int(M/2)-M_search), 1:N].shape)
        #print (indI)
    else:
        #Find the location of the maxima point in the lower part of the DHM hologram spectrum
        indI = np.unravel_index(np.argmax(pow_espec[1:M, 1:(int(N/2)-N_search)], axis=None), pow_espec[1:M, 1:(int(N/2)-N_search)].shape)
        #print (indI)
    
    #Determination of the ROI size. To do this, we use the theoretical size of the +1 or -1 diffraction orders according...
    #... to the distance of their centers to the DC coordinates (middle point in the DHM hologram spectrum).
    d = np.sqrt(np.power(indI[1]-M/2,2)+np.power(indI[0]-N/2,2))
    radius = d/3
    
    #Mask building
    mask = np.zeros((M,N))
    for j in range (M):
        for i in range (N):
            if np.power(j-indI[1], 2) + np.power(i-indI[0], 2) < np.power(radius, 2):
                mask[i,j] = 1
	
    tmp = field_spec*mask
    #field_spec.fill(0)
    
    #pow_espec = intensity(tmp, True)  
    #out = pow_espec
    
    #Coming back to spatial domain (retreiving filtered hologram)
    out = np.fft.ifftshift(tmp)
    out = np.fft.ifft2(out)
    out = np.fft.ifftshift(out)
	
    return out, indI, tmp