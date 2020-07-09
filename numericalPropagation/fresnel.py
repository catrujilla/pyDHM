
# libraries
import cv2
import numpy as np
from math import pi


def conversion_parameters(wave_length, distance):
    wave_length = (wave_length * 10 ** -9) / (1 * 10 ** -6)
    distance = distance / (1 * 10 ** -6)
    return wave_length, distance

# spatial filter
def spatial_filter(hologram, height, width):
    fft = np.fft.fft2(hologram)
    fft = np.fft.fftshift(fft)
    fft_holo_image = 20 * np.log(np.abs(fft))  # logaritm scale FFT
    minVal = np.amin(np.abs(fft_holo_image))
    maxVal = np.amax(np.abs(fft_holo_image))
    fft_holo_image = cv2.convertScaleAbs(fft_holo_image, alpha=255.0 / (maxVal - minVal),
                                         beta=-minVal * 255.0 / (maxVal - minVal))
    coordinates_ROI = cv2.selectROI("Fourier transform", fft_holo_image, fromCenter=True)  # module to  ROI
    x1_ROI = int(coordinates_ROI[1])
    y1_ROI = int(coordinates_ROI[0])
    x2_ROI = int(coordinates_ROI[1] + coordinates_ROI[3])
    y2_ROI = int(coordinates_ROI[0] + coordinates_ROI[2])
    ROI = fft[x1_ROI:x2_ROI, y1_ROI:y2_ROI]
    x = x2_ROI - x1_ROI
    y = y2_ROI - y1_ROI
    minX = int(width / 2 - x / 2)
    maxX = int(width / 2 + x / 2)
    minY = int(height / 2 - y / 2)
    maxY = int(height / 2 + y / 2)
    compensation = np.zeros((height, width), complex)
    compensation[minX:maxX, minY:maxY] = ROI
    compensation = np.fft.ifftshift(compensation)
    compensation = np.fft.ifft2(compensation)
    holo_filter = compensation

    #  holo_filter_image = cv2.magnitude(holo_filter[:, :, 0], holo_filter[:, :, 1])  #
    #  plt.imshow(holo_filter_image, cmap='gray'), plt.title('holo_filter')  # image in gray scale
    #  plt.show()  # show hologram
    return holo_filter


# Fresnel implementation
def fr(hologram, wave_length, distance, deltaX, deltaY):
    height, width = hologram.shape
    wave_length, distance = conversion_parameters(wave_length, distance)  # load parameters
    hologram = spatial_filter(hologram,height,width)
    # create a mesh-grid to operate in world-coordinates
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')

    # kernel to propagation
    kernel = np.exp2((-1j * pi / (wave_length * distance)) * (np.power(X * deltaX, 2) + np.power(Y * deltaY, 2)))

    # propagation
    reconstruction = hologram * kernel
    reconstruction = np.fft.fft2(reconstruction)
    reconstruction = np.fft.fftshift(reconstruction)
    reconstruction = (np.abs(reconstruction))
    return reconstruction

# Function to diffract a complex field using Fresnel approximation with
# Fourier method
# Inputs:
# field - complex field
# z - propagation distance
# wavelength - wavelength
# pitch_x/y - sampling pitches

def fresnel(field, z, wavelength, dx, dy):
    M, N = field.shape
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')

	dxout = (wavelength * z) / (M * dx)
	dyout = (wavelength * z) / (N * dy)
	
	k =  (2 * pi ) / wavelength
	
    z_phase = np.exp2((1j * k * z )/ (1i * wavelength * z))
	out_phase = np.exp2((1i * pi / (wavelength * z)) * (np.power(X * dxout, 2) + np.power(Y * dyout, 2)) )
	in_phase = np.exp2((1i * pi / (wavelength * z)) * (np.power(X * dx, 2) + np.power(Y * dy, 2)))

	tmp = (field * in_phase)
	tmp = np.fft.fftshift(tmp)
	tmp = np.fft.fft2(tmp)
	tmp = np.fft.fftshift(tmp)

	out = z_phase * out_phase * dx * dy * tmp
	
    return out

