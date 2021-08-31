# -*- coding: utf-8 -*-
"""
Title-->            PhaseCompensator definition script
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             21/02/2021
Last modified-->    31/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Script that implements the 3 phase compensators implemented in pyDHM
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

# libraries
import numpy as np
import math
from math import pi
import phaseShifting as phs

from matplotlib import pyplot as plt

def SOSR (inp0, inp1, inp2, inp3, upper, wavelength, deltaX, deltaY, s=1, steps=4):
    '''
    # Function to recover the phase information of a sample from 4 DHM slightly off-axis acquisitions (quadrature Phase Shifting method).
    # The idea was first proposed by Ferraro, we upgraded the method to accurately calculate the synthetic reference wave 
    # to obtain a fully compensated phase maps.
    # Inputs:
    # inpX - The input intensities (captured) pi/2 phase-shifted holograms
    # upper - Boolean variable defining whether the brightest point will be that on the upper part of the DHM hologram spectrum. 
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # deltaX, deltaY - Pixel dimensions of the camera sensor for the acquisition
    '''

    #Retrieving the input shape
    M, N = inp0.shape
    
    #Creating a mesh_grid to operate in world-coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY
    
    #Phase Shifting process
    
    #Location of the brightest point in the first phase-shifted DHM hologram spectrum
    #First, the need to Fourier transform and biuld the power spectrum
    field_spec = np.fft.fftshift(inp0)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)
    pow_espec = phs.intensity(field_spec, True)
    
    #Defining width and height of point of maxima in the ROI determination
    #10% is estimated to be the proper distance to separated an in-line hologram from an off-axis one in the Spectrum
    M_search = int(M-(M*0.1))
    N_search = int(N-(N*0.1))
    #print (M_search)
    
    #Determining the point of maxima in the power spectrum
    if (upper == True):
        #Find the location of the maxima point in the upper part of the DHM hologram spectrum
        indI = np.unravel_index(np.argmax(pow_espec[1:(int(M/2)-M_search), 1:N], axis=None), pow_espec[1:(int(M/2)-M_search), 1:N].shape)
        #print (indI)
    else:
        #Find the location of the maxima point in the lower part of the DHM hologram spectrum
        indI = np.unravel_index(np.argmax(pow_espec[1:M, 1:(int(N/2)-N_search)], axis=None), pow_espec[1:M, 1:(int(N/2)-N_search)].shape)
        #print (indI)
    
    #DC coordinates (centers of the DHM hologram Spectrums)
    fx_0, fy_0 = N/2, M/2
    
    #Coordinates of maxima points in the DHM hologram spectrums
    fx_1, fy_1 = indI[1]+1, indI[0]+1 # +1 is the correction of the offset added by numpy?
    
    #print (fx_0, fy_0, fx_1, fy_1)
    
    sum_max = 0
    arrayX = np.linspace(fx_1 - s, fx_1 + s, steps)
    arrayY = np.linspace(fy_1 - s, fy_1 + s, steps)

    k = (2 * pi) / wavelength
    
    print ("Synthetic building of the reference wave started....")
    
    for fx_temp in arrayX:
        for fy_temp in arrayY:
            #print (fx_temp, fy_temp)
            theta_x = math.asin((fx_0 - fx_temp) * wavelength / (N * deltaX))
            theta_y = math.asin((fy_0 - fy_temp) * wavelength / (M * deltaY))
            
            R0 = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
            R1 = 1j*np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
            R2 = (-1)*np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
            R3 = (-1j)*np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))

            # Four-Quarter Phase Shifting (Ferraro's proposal)
            Isynth = (R0*inp0) + (R1*inp1) + (R2*inp2) + (R3*inp3)
            
            #module 2-pi phase retrieval
            phase = np.angle(Isynth)

            #Thresholding process
            minVal = np.amin(phase)
            maxVal = np.amax(phase)
            phase_sca = (phase - minVal)/(maxVal - minVal)
            binary_phase = (phase_sca < 0.1)
            
            #plt.imshow(phase, cmap='gray'), plt.title("phase")  # image in gray scale
            #plt.show()  # show image
            
            #Applying the summation and thresholding metric
            sum = np.sum(np.sum(binary_phase))
            if sum > sum_max:
                x_max_out = fx_temp;
                y_max_out = fy_temp;
                sum_max = sum;
            #print (sum)
    
    #print(x_max_out, y_max_out)
    print ("Synthetic reference wave built.")
    
    #Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * deltaX))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * deltaY))
    R0 = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    R1 = 1j*np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    R2 = (-1)*np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    R3 = (-1j)*np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    
    # Four-Quarter Phase Shifting (Ferraro's proposal)
    Isynth = (R0*inp0) + (R1*inp1) + (R2*inp2) + (R3*inp3)
    
    comp_phase = Isynth
    
    print ("Phase compensation finished.")
    
    return comp_phase