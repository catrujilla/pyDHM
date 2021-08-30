# -*- coding: utf-8 -*-
"""
Title-->            Numerical propagator definition script
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             21/06/2021
Last modified-->    30/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Script that implements the 3 numerical propagator implemented in pyDHM
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

# libraries
import numpy as np
from math import pi

def fresnel(field, z, wavelength, dx, dy):
    """
    # Function to diffract a complex field using Fresnel approximation with
    # Fourier method
    # Inputs:
    # field - complex field
    # z - propagation distance
    # wavelength - wavelength
    # dx/dy - sampling pitches
    """
    M, N = field.shape
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')

    dxout = (wavelength * z) / (M * dx)
    dyout = (wavelength * z) / (N * dy)
    
    z_phase = np.exp2(1j * 2 * pi * z / wavelength) / (1j * wavelength * z)
   
    out_phase = np.exp2((1j * pi / (wavelength * z)) * (np.power(X * dxout, 2) + np.power(Y * dyout, 2)) )
    in_phase = np.exp2((1j * pi / (wavelength * z)) * (np.power(X * dx, 2) + np.power(Y * dy, 2)))

    tmp = (field * in_phase)
    tmp = np.fft.fftshift(tmp)
    tmp = np.fft.fft2(tmp)
    tmp = np.fft.fftshift(tmp)

    out = z_phase * out_phase * dx * dy * tmp
	
    return out

def bluestein(field, z, wavelength, dx, dy, dxout, dyout):
    """
    # Function to diffract a complex field using Fresnel approximation with 
    # Fresnel-Bluestein's method
    # Inputs:
    # field - complex field
    # z - propagation distance
    # lambda - wavelength
    # dx/dy - sampling pitches
    # dxout/dyout - output window pitches
    """
    
    M, N = field.shape
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')

    padx = int(M/2)
    pady = int(N/2)

    #print ('k: ' + str(k))
    #print ('wavelength: ' + str(wavelength))
    #print ('z: ' + str(z))
    
    z_phase = np.exp2(1j * 2 * pi * z / wavelength) / (1j * wavelength * z)
    
    #print ('z_phase: ' + str(z_phase))
    
    output_phase = np.exp2((-1j * pi/(wavelength*z)) * (dxout*(dx-dxout)*np.power(X, 2)+dyout*(dy-dyout)*np.power(Y, 2)))
    
    input_phase1 = np.exp2(((1j*pi)/(wavelength*z))*(dx*(dx-dxout)*np.power(X, 2)+dy*(dy-dyout)*np.power(Y, 2)))
    input_phase2 = np.exp2(((1j*pi)/(wavelength*z))*(dx*dxout*np.power(X, 2) + dy*dyout*np.power(Y, 2)))  
    
    f1 = np.pad(field*input_phase1, ((padx, padx), (pady, pady)), mode='constant')
    IP1 = np.fft.fftshift(f1)
    IP1 = np.fft.fft2(IP1)
    IP1 = np.fft.fftshift(IP1)
    
    f2 = np.pad(input_phase2, ((padx, padx), (pady, pady)), mode='constant')
    IP2 = np.fft.fftshift(f2)
    IP2 = np.fft.fft2(IP2)
    IP2 = np.fft.fftshift(IP2)  
    
    out = np.fft.ifftshift(IP1*IP2)
    out = np.fft.ifft2(out)
    out = np.fft.ifftshift(out)
    
    out = out[padx:padx+M,pady:pady+N]
    out = out*z_phase
    out = out*output_phase
    
    return out

def angularSpectrum(field, z, wavelength, dx, dy):
    '''
    # Function to diffract a complex field using the angular spectrum approach
    # Inputs:
    # field - complex field
    # z - propagation distance
    # wavelength - wavelength
    # dx/dy - sampling pitches
    '''
    M, N = field.shape
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')

    dfx = 1 / (dx * M)
    dfy = 1 / (dy * N)
    
    field_spec = np.fft.fftshift(field)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)
        
    phase = np.exp2(1j * z * pi * np.sqrt(np.power(1/wavelength, 2) - (np.power(X * dfx, 2) + np.power(Y * dfy, 2))))
	
    tmp = field_spec*phase
    
    out = np.fft.ifftshift(tmp)
    out = np.fft.ifft2(out)
    out = np.fft.ifftshift(out)
	
    return out