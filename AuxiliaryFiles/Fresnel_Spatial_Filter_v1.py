"""
Title-->            Fresnel
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             24/01/2019
Last modified-->    03/03/2020
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
Abstract -->        Algorithm that allow to elaborate the numerical reconstruction of digital holograms recorded in off
                    axis architecture trough Fresnel method and using an spatial filter to enhancement the
                    reconstruction.
"""

# libraries
import numpy as np
import cv2
from matplotlib import pyplot as plt
from math import pi
from PIL import Image


# main fuction
def main():
    holo = read_hologram()  # read the hologram (input_Field)
    height, width = holo.shape  # get size of the image
    wave_length, distance, dimensionX, dimensionY, deltaX, deltaY = load_parameters(height, width)  # load parameters
    holo_filter = spatial_filter(holo, height, width)
    reconstruction = fresnel(holo_filter, height, width, wave_length, distance, deltaX, deltaY)
    intensity(reconstruction)


# read hologram
def read_hologram():
    namefile = 'mask.jpg'
    #namefile = 'die_1.jpg'
    holo = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
    plt.imshow(holo, cmap='gray'), plt.title('hologram')  # image in gray scale
    plt.show()  # show hologram
    return holo


# reconstruction parameters (units -- > micrometers)
def load_parameters(height, width):
    wave_length = (632.8 * 10 ** -9) / (1 * 10 ** -6)  # wavelength used to the digital reconstruction
    distance = 1.05 / (1 * 10 ** -6)  # distance of propagation
    dimensionX = (5.2 * 10 ** -3) / (1 * 10 ** -6)  # X dimension of the sensor
    dimensionY = (5.2 * 10 ** -3) / (1 * 10 ** -6)  # Y dimension of the sensor
    deltaX = dimensionX / height  # pixel X dimension
    deltaY = dimensionY / width  # pixel Y dimension
    return wave_length, distance, dimensionX, dimensionY, deltaX, deltaY


# spatial filter
def spatial_filter(holo, height, width):
    holoFT = np.float32(holo)  # convertion of data to float
    fft_holo = cv2.dft(holoFT, flags=cv2.DFT_COMPLEX_OUTPUT)  # FFT of hologram
    fft_holo = np.fft.fftshift(fft_holo)
    fft_holo_image = 20 * np.log(cv2.magnitude(fft_holo[:, :, 0], fft_holo[:, :, 1]))  # logaritm scale FFT
    minVal = np.amin(np.abs(fft_holo_image))
    maxVal = np.amax(np.abs(fft_holo_image))
    fft_holo_image = cv2.convertScaleAbs(fft_holo_image, alpha=255.0 / (maxVal - minVal),
                                         beta=-minVal * 255.0 / (maxVal - minVal))
    coordinates_ROI = cv2.selectROI("Fourier transform", fft_holo_image, fromCenter=True)  # module to  ROI
    x1_ROI = int(coordinates_ROI[1])
    y1_ROI = int(coordinates_ROI[0])
    x2_ROI = int(coordinates_ROI[1] + coordinates_ROI[3])
    y2_ROI = int(coordinates_ROI[0] + coordinates_ROI[2])
    ROI = fft_holo[x1_ROI:x2_ROI, y1_ROI:y2_ROI]
    x = x2_ROI - x1_ROI
    y = y2_ROI - y1_ROI
    minX = int(width / 2 - x / 2)
    maxX = int(width / 2 + x / 2)
    minY = int(height / 2 - y / 2)
    maxY = int(height / 2 + y / 2)
    compensation = np.zeros((height, width, 2))
    compensation[minX:maxX, minY:maxY] = ROI
    #  compensation_image = cv2.magnitude(compensation[:, :, 0], compensation[:, :, 1])  #
    #  plt.imshow(compensation_image, cmap='gray'), plt.title('Spatial Frecuencies of interest')  # image in gray scale
    #  plt.show()  # show spatial filter
    compensation = np.fft.ifftshift(compensation)
    holo_filter = cv2.idft(compensation, flags=cv2.DFT_INVERSE)
    #  holo_filter_image = cv2.magnitude(holo_filter[:, :, 0], holo_filter[:, :, 1])  #
    #  plt.imshow(holo_filter_image, cmap='gray'), plt.title('holo_filter')  # image in gray scale
    #  plt.show()  # show hologram
    return holo_filter


# Fresnel reconstruction
def fresnel(holo_filter, height, width, wave_length, distance, deltaX, deltaY):
    # create a mesh-grid to operate in world-coordinates
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY
    # kernel to propagation
    kernel = np.exp2((-1j * pi / (wave_length * distance)) * (np.power(X * deltaX, 2) + np.power(Y * deltaY, 2)))
    #  propagation
    holo_filter_real = holo_filter[:, :, 0]
    holo_filter_imag = holo_filter[:, :, 1]
    holo_filter = np.zeros((height, width), complex)
    for p in range(height):
        for q in range(width):
            holo_filter[p, q] = complex(holo_filter_real[p, q], holo_filter_imag[p, q])
    reconstruction = holo_filter * kernel
    reconstruction_real = reconstruction.real
    reconstruction_imag = reconstruction.imag
    reconstruction = np.dstack([reconstruction_real, reconstruction_imag])
    reconstruction = cv2.dft(reconstruction, flags=cv2.DFT_COMPLEX_OUTPUT)
    reconstruction = np.fft.fftshift(reconstruction)
    return reconstruction


def intensity(reconstruction):
    modu_reconstruction = (cv2.magnitude(reconstruction[:, :, 0], reconstruction[:, :, 1]))
    plt.imshow(modu_reconstruction, cmap='gray'), plt.title('Reconstruction')  # image in gray scale
    plt.show()  # show hologram
    saveImage(modu_reconstruction, 'Reconstruction.jpg')


def saveImage(image, title):
    minVal = np.amin(image)
    maxVal = np.amax(image)
    image = cv2.convertScaleAbs(image, alpha=255.0 / (maxVal - minVal),
                                beta=-minVal * 255.0 / (maxVal - minVal))
    image = Image.fromarray(image)
    image.save(title)


main()
