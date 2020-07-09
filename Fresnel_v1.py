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
#import cv2
from matplotlib import pyplot as plt
from math import pi
from PIL import Image

# main fuction
def main():
    namefile = 'mask.jpg'
    holo = read_hologram(namefile)  # read the hologram (input_Field)
    height, width = holo.shape  # get size of the image
    wave_length, distance, dimensionX, dimensionY, deltaX, deltaY = load_parameters(height, width)  # load parameters
    reconstruction = fresnel(holo, height, width, wave_length, distance, deltaX, deltaY)
    intensity(reconstruction)


# read hologram
def read_hologram(namefile):
    #holo = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
    holo = Image.open(namefile) 
    plt.imshow(holo, cmap='gray'), plt.title('hologram')  # image in gray scale
    plt.show()  # show hologram
    return holo


# reconstruction parameters (units -- > micrometers)
def load_parameters(height, width):
    wave_length = (632.8*10**-9)/(1*10**-6)  # wavelength used to the digital reconstruction
    distance = 0.05/(1*10**-6)  # distance of propagation
    dimensionX = (5.2*10**-3)/(1*10**-6)  # X dimension of the sensor
    dimensionY = (5.2*10**-3)/(1*10**-6)  # Y dimension of the sensor
    deltaX = dimensionX / height  # pixel X dimension
    deltaY = dimensionY / width  # pixel Y dimension
    return wave_length, distance, dimensionX, dimensionY, deltaX, deltaY


# Fresnel reconstruction
def fresnel(holo, height, width, wave_length, distance, deltaX, deltaY):
    # create a mesh-grid to operate in world-coordinates
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY
    # kernel to propagation
    kernel = np.exp2((-1j * pi / (wave_length * distance)) * (np.power(X * deltaX, 2) + np.power(Y * deltaY, 2)))
    #propagation
    reconstruction = holo * kernel
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
