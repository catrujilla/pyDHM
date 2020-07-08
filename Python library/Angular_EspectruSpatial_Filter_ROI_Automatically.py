"""
Title-->            Angular Espectrum - Image Plane holograms
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             03/03/2019
Last modified-->    24/01/2020
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
Abstract -->        Algorithm that allow to elaborate the numerical reconstruction of digital holograms recorded in
                    image plane architecture and off-axis setup, trough spatial Filter method
Links-->          - https://medium.com/@y1017c121y/python-computer-vision-tutorials-image-fourier-transform-part-2-ec980
                    3e63993 (FFT time)
"""

# libraries
import numpy as np
import cv2
import math
from matplotlib import pyplot as plt
from math import pi


def main():
    holo = read_hologram()
    height, width = holo.shape  # get size of the image
    wave_length, distance, dimensionX, dimensionY, deltaX, deltaY = load_parameters(height, width)
    holo_filter, indI = spatial_filter(holo, height, width)
    x_max_out, y_max_out = automatically_angular_spectrum(holo_filter, indI, height, width, wave_length, deltaX, deltaY)
    reconstruction = best_reconstruction(holo_filter, x_max_out, y_max_out, height, width, wave_length, deltaX, deltaY)
    amplitude(reconstruction)
    phase(reconstruction, height, width)


# read hologram
def read_hologram():
    namefile = 'USAF_20x.jpg'
    holo = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
    imageShow(holo, 'hologram')
    return holo


# reconstruction parameters (units -- > micro-meters)
def load_parameters(height, width):
    wave_length = 0.633  # wavelength used to the digital reconstruction
    distance = 0  # distance of propagation
    dimensionX = 7065.6  # X dimension of the sensor
    dimensionY = 7065.6  # Y dimension of the sensor
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
    imageShow(fft_holo_image, 'FFT hologram')
    indI = np.unravel_index(np.argmax(fft_holo_image[1:500, 1:500], axis=None), fft_holo_image[1:500, 1:500].shape)
    coordinates_ROI = cv2.selectROI(fft_holo_image, fromCenter=True)  # module to  ROI
    imCrop = fft_holo_image[int(coordinates_ROI[1]):int(coordinates_ROI[1] + coordinates_ROI[3]),
             int(coordinates_ROI[0]):int(coordinates_ROI[0] + coordinates_ROI[2])]
    x1_ROI = int(coordinates_ROI[1])
    y1_ROI = int(coordinates_ROI[0])
    x2_ROI = int(coordinates_ROI[1] + coordinates_ROI[3])
    y2_ROI = int(coordinates_ROI[0] + coordinates_ROI[2])
    holo_filter = np.zeros((height, width, 2))
    holo_filter[x1_ROI:x2_ROI, y1_ROI: y2_ROI] = 1
    holo_filter = holo_filter * fft_holo
    holo_filter = np.fft.ifftshift(holo_filter)
    holo_filter = cv2.idft(holo_filter, flags=cv2.DFT_INVERSE)
    filterMatrix = (cv2.magnitude(holo_filter[:, :, 0], holo_filter[:, :, 1]))  # logaritm scale for the FFT
    minVal = np.amin(np.abs(filterMatrix))
    maxVal = np.amax(np.abs(filterMatrix))
    filterMatrix = cv2.convertScaleAbs(filterMatrix, alpha=255.0 / (maxVal - minVal),
                                       beta=-minVal * 255.0 / (maxVal - minVal))
    imageShow(filterMatrix, 'Holo Filter')
    return holo_filter, indI


#  Image Plane reconstruction
def automatically_angular_spectrum(holo_filter, indI, height, width, wave_length, deltaX, deltaY):
    # create a mesh_grid to operate in world-coordinates
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY

    holo_filter_real = holo_filter[:, :, 0]
    holo_filter_imag = holo_filter[:, :, 1]
    holo_filter = np.zeros((height, width), complex)
    for p in range(height):
        for q in range(width):
            holo_filter[p, q] = complex(holo_filter_real[p, q], holo_filter_imag[p, q])
    #  kernel to propagation
    fx_0 = height / 2
    fy_0 = width / 2
    #fx_1 = 268.5000
    #fy_1 = 287.4000
    fx_1 = indI[1]
    fy_1 = indI[0]

    s = 1
    step = 0.5
    i = 0
    sum_max = 0
    arrayX = np.linspace(fx_1 - s, fx_1 + s, 21)
    arrayY = np.linspace(fy_1 - s, fy_1 + s, 21)
    #array = [268.5000, 287.4000]
    #  for fx_temp
    k = (2 * pi) / wave_length
    for fx_temp in arrayX:
        for fy_temp in arrayY:
            theta_x = math.asin((fx_0 - fx_temp) * wave_length / (height * deltaX))
            theta_y = math.asin((fy_0 - fy_temp) * wave_length / (width * deltaY))
            ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
            #  propagation
            reconstruction = holo_filter * ref_wave
            reconstruction_real = reconstruction.real
            reconstruction_imag = reconstruction.imag
            reconstruction = np.dstack([reconstruction_real, reconstruction_imag])
            reconstruction_real = reconstruction[:, :, 0]
            reconstruction_imag = reconstruction[:, :, 1]
            reconstruction = np.zeros((height, width), complex)
            for p in range(height):
                for q in range(width):
                    reconstruction[p, q] = complex(reconstruction_real[p, q], reconstruction_imag[p, q])
            phase = np.angle(reconstruction)
            minVal = np.amin(phase)
            maxVal = np.amax(phase)
            phase = cv2.convertScaleAbs(phase, alpha=255.0 / (maxVal - minVal),
                                        beta=-minVal * 255.0 / (maxVal - minVal))
            binary_phase = cv2.threshold(phase, 50, 255, cv2.THRESH_BINARY)
            sum = np.sum(np.sum(binary_phase))
            if sum > sum_max:
                x_max_out = fx_temp;
                y_max_out = fy_temp;
                sum_max = sum;
            #print(sum)
            #plt.imshow(phase, 'gray'), plt.title('phase')  # image in gray scale
            #plt.show()
    print(x_max_out, y_max_out)
    return x_max_out, y_max_out

def best_reconstruction(holo_filter, x_max_out, y_max_out, height, width, wave_length, deltaX, deltaY):
    # create a mesh_grid to operate in world-coordinates
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY
    fx_0 = height / 2
    fy_0 = width / 2
    k = (2 * pi) / (wave_length)
    theta_x = math.asin((fx_0 - x_max_out) * wave_length / (height * deltaX))
    theta_y = math.asin((fy_0 - y_max_out) * wave_length / (width * deltaY))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    #  propagation
    holo_filter_real = holo_filter[:, :, 0]
    holo_filter_imag = holo_filter[:, :, 1]
    holo_filter = np.zeros((height, width), complex)
    for p in range(height):
        for q in range(width):
            holo_filter[p, q] = complex(holo_filter_real[p, q], holo_filter_imag[p, q])
    reconstruction = holo_filter * ref_wave
    reconstruction_real = reconstruction.real
    reconstruction_imag = reconstruction.imag
    reconstruction = np.dstack([reconstruction_real, reconstruction_imag])
    return reconstruction

def amplitude(reconstruction):
    amplitude = (cv2.magnitude(reconstruction[:, :, 0], reconstruction[:, :, 1]))  # logaritm scale for the FFT
    plt.imshow(amplitude, cmap='gray'), plt.title('Amplitude')  # image in gray scale
    plt.show()


def phase(reconstruction, height, width):
    reconstruction_real = reconstruction[:, :, 0]
    reconstruction_imag = reconstruction[:, :, 1]
    reconstruction = np.zeros((height, width), complex)
    for p in range(height):
        for q in range(width):
            reconstruction[p, q] = complex(reconstruction_real[p, q], reconstruction_imag[p, q])
    phase = np.angle(reconstruction)
    plt.imshow(phase, cmap='gray'), plt.title('phase')  # image in gray scale
    plt.show()


def imageShow(image, title):
    plt.imshow(image, cmap='gray'), plt.title(title)  # image in gray scale
    plt.show()  # show hologram


main()
