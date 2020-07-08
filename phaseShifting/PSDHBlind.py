"""
Title-->            PS-DHM blind method
Author-->           Raul Castaneda*, Ana Doblas*, Carlos Trujillo**
                    *University of Memphis, USA; Optical Imaging Research lab (OIRL)
                    **EAFIT, Colombia; 
Date-->             04/06/2020
Last modified-->    04/06/2020

Abstract -->        Library to  Algorithm that allow reconstruction of numerical reconstruction for digital holograms
                    recorded in
                    architecture, image plane, off-axis setup and in telecentric configuration
"""

import numpy as np
import math
from scipy.optimize import minimize


# cut the holograms into a square shape
def cut_holograms2(hologram1, hologram2):
    height, width = hologram1.shape
    if width > height:
        cut = round(int(width - height) / 2)
        hologram1 = hologram1[0:height, cut:(width - cut)]
        hologram2 = hologram2[0:height, cut:(width - cut)]
    elif width < height:
        cut = round(int(height - width) / 2)
        hologram1 = hologram1[cut:(height - cut), 0:width]
        hologram2 = hologram2[cut:(height - cut), 0:width]
    else:
        hologram1 = hologram1
        hologram2 = hologram2
    return hologram1, hologram2


# cut the holograms into a square shape
def cut_holograms3(hologram1, hologram2, hologram3):
    height, width = hologram1.shape
    if width > height:
        cut = round(int(width - height) / 2)
        hologram1 = hologram1[0:height, cut:(width - cut)]
        hologram2 = hologram2[0:height, cut:(width - cut)]
        hologram3 = hologram3[0:height, cut:(width - cut)]
    elif width < height:
        cut = round(int(height - width) / 2)
        hologram1 = hologram1[cut:(height - cut), 0:width]
        hologram2 = hologram2[cut:(height - cut), 0:width]
        hologram3 = hologram3[cut:(height - cut), 0:width]
    else:
        hologram1 = hologram1
        hologram2 = hologram2
        hologram3 = hologram3
    return hologram1, hologram2, hologram3


# find the pixel position for the max values
def find_peaks(hologram1, height, width):
    # FFT of the hologram
    fft = np.fft.fft2(hologram1)
    fft = np.fft.fftshift(fft)

    # Find the max peaks I and II quadrants
    height_half = round(int(height / 2))
    maskII = np.zeros((height, width))
    maskII[0:height_half - 1, 0:width] = 1
    fft_images_II = fft * maskII
    maximunII = np.amax(fft_images_II)
    fy_max_d1, fx_max_d1 = np.where(fft_images_II == maximunII)

    # Find the max peaks III and IV quadrants
    maskIV = np.zeros((height, width))
    maskIV[height_half + 1:height, 0:width] = 1
    fft_images_IV = fft * maskIV
    maximun2 = np.amax(fft_images_IV)
    fy_max_d2, fx_max_d2 = np.where(fft_images_IV == maximun2)
    return fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2


# cost function for two holograms
def cost_function2(test_tetha, hologram1, hologram2, fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2):
    theta1 = test_tetha[0]
    theta2 = test_tetha[1]
    ones = np.ones((2, 2))
    M = [[1, np.exp(1j*theta1)],
         [1, np.exp(1j*theta2)]]
    M = M * ones
    Minv = np.linalg.inv(M)
    d1 = Minv[0, 0]*hologram1 + Minv[0, 1]*hologram2

    # FFT d1
    FTd1 = np.fft.fft2(d1)
    FTd1 = np.fft.fftshift(FTd1)

    # calculate J
    J = 1 - (np.abs(FTd1[fx_max_d1, fy_max_d1]) - np.abs(FTd1[fx_max_d2, fy_max_d2])) / \
        (np.abs(FTd1[fx_max_d2, fy_max_d2]) + np.abs(FTd1[fx_max_d1, fy_max_d1]))
    J = J[0]
    return J


# cost function for two holograms
def cost_function3(test_tetha, hologram1, hologram2, hologram3, fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2):
    theta1 = test_tetha[0]
    theta2 = test_tetha[1]
    theta3 = test_tetha[2]
    ones = np.ones((3, 3))
    M = [[1, np.exp(1j*theta1), np.exp(-1j*theta1)],
         [1, np.exp(1j*theta2), np.exp(-1j*theta2)],
         [1, np.exp(1j*theta3), np.exp(-1j*theta3)]]
    M = M * ones
    Minv = np.linalg.inv(M)
    d3 = Minv[2, 0]*hologram1 + Minv[2, 1]*hologram2 + Minv[2, 2]*hologram3

    # FFT d3
    FTd3 = np.fft.fft2(d3)
    FTd3 = np.fft.fftshift(FTd3)

    # calculate J
    J = 1 + (1 - (np.abs(FTd3[fx_max_d2, fy_max_d2]) - np.abs(FTd3[fx_max_d1, fy_max_d1])) /
             (np.abs(FTd3[fx_max_d2, fy_max_d2]) + np.abs(FTd3[fx_max_d1, fy_max_d1])))
    J = J[0]
    return J


# calculate d1 and d2 for two holograms
def calculate_d1_d2(x, hologram1, hologram2):
    theta1 = x[0]
    theta2 = x[1]
    ones = np.ones((2, 2))
    M = [[1, np.exp(1j*theta1)],
         [1, np.exp(1j*theta2)]]
    M = M * ones
    Minv = np.linalg.inv(M)
    d2 = Minv[1, 0]*hologram1 + Minv[1, 1]*hologram2

    # FFT d2
    FTd2 = np.fft.fft2(d2)
    FTd2 = np.fft.fftshift(FTd2)
    return FTd2, d2


# calculate d1, d2 and d3 for three holograms
def calculate_d1_d2_d3(x, hologram1, hologram2, hologram3):
    theta1 = x[0]
    theta2 = x[1]
    theta3 = x[2]
    ones = np.ones((3, 3))
    M = [[1, np.exp(1j*theta1), np.exp(-1j*theta1)],
         [1, np.exp(1j*theta2), np.exp(-1j*theta2)],
         [1, np.exp(1j*theta3), np.exp(-1j*theta3)]]
    M = M * ones
    Minv = np.linalg.inv(M)
    d1 = Minv[0, 0]*hologram1 + Minv[0, 1]*hologram2 + Minv[0, 2]*hologram3
    d2 = Minv[1, 0]*hologram1 + Minv[1, 1]*hologram2 + Minv[1, 2]*hologram3
    d3 = Minv[2, 0]*hologram1 + Minv[2, 1]*hologram2 + Minv[2, 2]*hologram3

    # FFT d2
    FTd2 = np.fft.fft2(d2)
    FTd2 = np.fft.fftshift(FTd2)
    return FTd2, d2


# spatial filter for two holograms
def spatial_filter_d2(height, width, FTd2):
    # Find the max peaks I or II quadrant
    height_half = round(int(height / 2))
    mask = np.zeros((height, width))
    mask[0:height_half - 1, 0:width] = 1
    fft_images_II = FTd2 * mask
    maximum = np.amax(fft_images_II)
    fy_max, fx_max = np.where(fft_images_II == maximum)
    # image_show(fft_images_II, 'II quadrant'), print(fx_max, fy_max)

    # Filter hologram
    '''
    radius = 100
    for p in range(0, height):
        for q in range(0, width):
            if np.sqrt((p - fy_max) ** 2 + (q - fx_max) ** 2) < radius:
                mask[p, q] = 1
    '''
    limit = 100
    x1 = round(int(fx_max - limit))
    x2 = round(int(fx_max + limit))
    y1 = round(int(fy_max - limit))
    y2 = round(int(fy_max + limit))
    mask[x1:x2, y1:y2] = 1

    hologram_filter = FTd2 * mask
    hologram_filter = np.fft.ifftshift(hologram_filter)
    hologram_filter = np.fft.ifft2(hologram_filter)
    return hologram_filter


# calculate the reference wave
def reference_wave(height, width, wave_length, deltaX, deltaY, fx_max, fy_max):
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY
    fx_0 = height / 2
    fy_0 = width / 2
    k = (2 * math.pi) / wave_length
    theta_x = math.asin((fx_0 - fx_max) * wave_length / (height * deltaX))
    theta_y = math.asin((fy_0 - fy_max) * wave_length / (width * deltaY))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    return ref_wave


# numerical phase reconstruction
def numerical_reconstructions(d2, ref_wave):
    reconstruction = d2 * ref_wave
    phase = np.angle(reconstruction)
    return phase


def raw2(hologram1, hologram2, wave_length, deltaX, deltaY):
    hologram1, hologram2 = cut_holograms2(hologram1, hologram2)
    height, width = hologram1.shape
    fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2 = find_peaks(hologram1, height, width)
    test_tetha = np.random.randint(0, 360, (1, 2))
    test_tetha = test_tetha * math.pi / 180
    test_tetha = [test_tetha[0, 0], test_tetha[0, 1]]
    res = minimize(cost_function2, test_tetha, args=(hologram1, hologram2, fx_max_d1, fy_max_d1, fx_max_d2,
                                                    fy_max_d2), method='Cobyla', tol=1e-6)
    x = res.x
    FTd2, d2 = calculate_d1_d2(x, hologram1, hologram2)
    hologram_filter = spatial_filter_d2(height, width, FTd2)
    ref_wave = reference_wave(height, width, wave_length, deltaX, deltaY, fx_max_d1, fy_max_d1)
    phase = numerical_reconstructions(hologram_filter, ref_wave)
    return phase


def raw3(hologram1, hologram2, hologram3,  wave_length, deltaX, deltaY):
    hologram1, hologram2, hologram3 = cut_holograms3(hologram1, hologram2, hologram3)
    height, width = hologram1.shape
    fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2 = find_peaks(hologram1, height, width)
    test_tetha = np.random.randint(0, 360, (1, 3))
    test_tetha = test_tetha * math.pi / 180
    test_tetha = [test_tetha[0, 0], test_tetha[0, 1], test_tetha[0, 2]]
    res = minimize(cost_function3, test_tetha, args=(hologram1, hologram2, hologram3, fx_max_d1, fy_max_d1, fx_max_d2,
                                                    fy_max_d2), method='Cobyla', tol=1e-6)
    x = res.x
    FTd2, d2 = calculate_d1_d2_d3(x, hologram1, hologram2, hologram3)
    #hologram_filter = spatial_filter_d2(height, width, FTd2)
    ref_wave = reference_wave(height, width, wave_length, deltaX, deltaY, fx_max_d1, fy_max_d1)
    phase = numerical_reconstructions(d2, ref_wave)
    return phase

