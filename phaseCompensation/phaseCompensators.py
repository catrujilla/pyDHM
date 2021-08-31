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
Links-->          - https://github.com/catrujilla/pyDHM
"""

# libraries
import numpy as np
import math
from math import pi
import phaseCompensation as phc

from matplotlib import pyplot as plt


def FRS (inp, upper, wavelength, deltaX, deltaY, s=2, steps=10):
    '''
    # Function to compensate phase maps of off-axis DHM via the full ROI search algorithm.
    # Inputs:
    # inp - The input intensity (captured) hologram
    # upper - Boolean variable defining whether the ROI will be that on the upper part of the DHM hologram spectrum. 
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # deltaX, deltaY - Pixel dimensions of the camera sensor for the acquisition
    '''

    #Retrieving the input shape
    M, N = inp.shape
    
    #Creating a mesh_grid to operate in world-coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY
    
    #The spatial filtering process is executed
    print ("Spatial filtering process started.....")
    holo_filter, indI, tmp = phc.spatialFiltering(inp, M, N, upper)
    print ("Spatial filtering process finished.")
    
    #out = phc.intensity(tmp, True)
    #plt.imshow(out, cmap='gray'), plt.title("holo_filterFreq")  # image in gray scale
    #plt.show()  # show image
    
    #Phase Compensation process
    
    #DC coordinates (center of DHM hologram Spectrum)
    fx_0, fy_0 = N/2, M/2
    
    #Coordintes of maxima in the DHM hologram spectrum
    fx_1, fy_1 = indI[1]+1, indI[0]+1 # +1 is the correction of the offset added by numpy?
    
    #print (fx_0, fy_0, fx_1, fy_1)
    
    sum_max = 0
    arrayX = np.linspace(fx_1 - s, fx_1 + s, steps)
    arrayY = np.linspace(fy_1 - s, fy_1 + s, steps)

    k = (2 * pi) / wavelength
    
    print ("Phase compensation started....")
    
    for fx_temp in arrayX:
        for fy_temp in arrayY:
            #print (fx_temp, fy_temp)
            theta_x = math.asin((fx_0 - fx_temp) * wavelength / (N * deltaX))
            theta_y = math.asin((fy_0 - fy_temp) * wavelength / (M * deltaY))
            ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
            #Compensation of the tilting angle for the off-axis acquisition
            reconstruction = holo_filter * ref_wave
            #module 2-pi phase retrieval
            phase = np.angle(reconstruction)
            #Thresholding process
            minVal = np.amin(phase)
            maxVal = np.amax(phase)
            phase_sca = (phase - minVal)/(maxVal - minVal)
            binary_phase = (phase_sca > 0.2)
            
            #plt.imshow(binary_phase, cmap='gray'), plt.title("binary_phase")  # image in gray scale
            #plt.show()  # show image
            
            #Applying the summation and thresholding metric
            sum = np.sum(np.sum(binary_phase))
            if sum > sum_max:
                x_max_out = fx_temp;
                y_max_out = fy_temp;
                sum_max = sum;
            #print (sum)

    #Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * deltaX))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * deltaY))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    comp_phase = holo_filter * ref_wave
    
    print ("Phase compensation finished.")
    #print(x_max_out, y_max_out)
    
    return comp_phase
    
def FastRS (inp, upper, wavelength, deltaX, deltaY, s=5, step=0.2):
    '''
    # Function to compensate phase maps of off-axis DHM via the fast ROI search algorithm.
    # Inputs:
    # inp - The input intensity (captured) hologram
    # upper - Boolean variable defining whether the ROI will be that on the upper part of the DHM hologram spectrum. 
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # deltaX, deltaY - Pixel dimensions of the camera sensor for the acquisition
    '''

    #Retrieving the input shape
    M, N = inp.shape
    
    #Creating a mesh_grid to operate in world-coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY
    
    #The spatial filtering process is executed
    print ("Spatial filtering process started.....")
    holo_filter, indI, tmp = phc.spatialFiltering(inp, M, N, upper)
    print ("Spatial filtering process finished.")
    
    #out = phc.intensity(tmp, True)
    #plt.imshow(out, cmap='gray'), plt.title("holo_filterFreq")  # image in gray scale
    #plt.show()  # show image
    
    #Phase Compensation process
    
    #DC coordinates (center of DHM hologram Spectrum)
    fx_0, fy_0 = N/2, M/2
    
    #Coordintes of maxima in the DHM hologram spectrum
    fx_1, fy_1 = indI[1]+1, indI[0]+1 # +1 is the correction of the offset added by numpy?
    
    #print (fx_0, fy_0, fx_1, fy_1)
    
    print ("Phase compensation started....")

    k = (2 * pi) / wavelength

    fin = 0
    sum_max = 0

    fx = fx_1
    fy = fy_1

    G_temp = s
        
    while fin == 0:
  
        sum_max=0 #small number for the metric (thresholding)    
    
        arrayY = np.linspace( int(10*(fy-step*G_temp)), int(10*(fy+step*G_temp)) , int(10*step))
        arrayX = np.linspace( int(10*(fx-step*G_temp)), int(10*(fx+step*G_temp)) , int(10*step))
    
        for fx_temp in arrayX:
            for fy_temp in arrayY:
                
                fx_tmp, fy_tmp = fx_temp/10, fy_temp/10
                
                #print (fx_tmp, fy_tmp)
                
                theta_x = math.asin((fx_0 - fx_tmp) * wavelength / (N * deltaX))
                theta_y = math.asin((fy_0 - fy_tmp) * wavelength / (M * deltaY))
                ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))

                #Compensation of the tilting angle for the off-axis acquisition
                reconstruction = holo_filter * ref_wave
                #module 2-pi phase retrieval
                phase = np.angle(reconstruction)
                #Thresholding process
                minVal = np.amin(phase)
                maxVal = np.amax(phase)
                phase_sca = (phase - minVal)/(maxVal - minVal)
                binary_phase = (phase_sca > 0.2)
            
                #plt.imshow(phase, cmap='gray'), plt.title("binary_phase")  # image in gray scale
                #plt.show()  # show image
            
                #Applying the summation and thresholding metric
                sum = np.sum(np.sum(binary_phase))
                if (sum > sum_max):
                    x_max_out = fx_tmp
                    y_max_out = fy_tmp
                    sum_max=sum
                #print (sum)
  
        G_temp = G_temp - 1
  
        if x_max_out == fx and y_max_out == fy:
            fin = 1;       
  
        fx = x_max_out
        fy = y_max_out
        
    #print(x_max_out, y_max_out)

    #Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * deltaX))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * deltaY))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    comp_phase = holo_filter * ref_wave
    
    print ("Phase compensation finished.")
    
    return comp_phase