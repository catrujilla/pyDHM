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
	
    k =  (2 * pi ) / wavelength
	
    z_phase = np.exp2((1j * k * z )/ (1j * wavelength * z))
    out_phase = np.exp2((1j * pi / (wavelength * z)) * (np.power(X * dxout, 2) + np.power(Y * dyout, 2)) )
    in_phase = np.exp2((1j * pi / (wavelength * z)) * (np.power(X * dx, 2) + np.power(Y * dy, 2)))

    tmp = (field * in_phase)
    tmp = np.fft.fftshift(tmp)
    tmp = np.fft.fft2(tmp)
    tmp = np.fft.fftshift(tmp)

    out = z_phase * out_phase * dx * dy * tmp
	
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