# -*- coding: utf-8 -*-
"""
Title-->            PhaseCompensator definition script
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             21/06/2021
Last modified-->    30/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Script that implements the 3 phase compensators implemented in pyDHM
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

# libraries
import numpy as np
#from math import pi


def FRS (inp):
    '''
    # Function to spatial-filter a DHM hologram prior reconstruction
    # Inputs:
    # inp - The input intensity (captured) hologram
    '''

    M, N = inp.shape
    field_spec = np.fft.fftshift(inp)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)

    out = intensity (field_spec, True)  

    imageShow (out, "Espectro")
        
    #Location and size of the ROI for the +a diffraction order info
    x0 = 235
    y0 = 650
    radius = 100
    
    #Mask building
    mask = np.zeros((M,N))
    for j in range (M):
        for i in range (N):
            if np.power(j-x0, 2) + np.power(i-y0, 2) < np.power(radius, 2):
                mask[i,j] = 1
	
    tmp = field_spec*mask
    field_spec.fill(0)

    #Centering the +1 diffraction order information
    for j in range (M):
        for i in range (N):
            new_i = i-y0+int(N/2)
            new_j = j-x0+int(M/2)
            if (new_i<960 and new_j<960 and new_i>0 and new_j>0):
                field_spec[new_i,new_j] = tmp[i,j]
    
    out = np.fft.ifftshift(field_spec)
    out = np.fft.ifft2(out)
    out = np.fft.ifftshift(out)
        
    return out