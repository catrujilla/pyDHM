# -*- coding: utf-8 -*-
"""
Title-->            Display options utility script
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             03/03/2019
Last modified-->    16/07/2020
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Script that implements the different methods to render the resulting complex field data
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

import numpy as np
from matplotlib import pyplot as plt

# Function to calcule the amplitude representation of a given complex field
# Inputs:
# inp - The input complex field
# log - boolean variable to determine if a log representation is applied
def amplitude(inp, log):
    out = np.abs(inp)

    if log == True:
        out = 20 * np.log(out)

    return out


# Function to calcule the intensity representation of a given complex field
# Inputs:
# inp - The input complex field
# log - boolean variable to determine if a log representation is applied
def intensity(inp, log):
    out = np.abs(inp)
    out = out * out

    if log == True:
        out = 20 * np.log(out)
        out[out == np.inf] = 0
        out[out == -np.inf] = 0

    return out


# Function to calcule the phase representation of a given complex field using the
# function 'angle'
# Inputs:
# inp - The input complex field
def phase(inp):
    out = np.angle(inp)

    return out


# Function to calcule the Fourier transform of a given field using the
# 'fft of numpy'
# Inputs:
# inp - The input field
def FT(inp):
    FT = np.fft.fft2(inp)
    FT = np.fft.fftshift(FT)



    return FT


def sfc(M, N, field, radius, centX, centY):
    """
    # Function to create a spatial filter with a circle mask
    # Inputs:
    # M - Pixel size x
    # N - Pixel size y
    # field - The field to be filtered
    # radius - dimension of the radius for the circle
    # centX - center circle position in x axis
    # centY - center circle position in y axis
    """
    circle = np.zeros((radius * 2, radius * 2), dtype=int)
    for p in range(0, radius * 2):
        for q in range(0, radius * 2):
            if np.sqrt((p - radius) ** 2 + (q - radius) ** 2) < radius:
                circle[p, q] = 1

    (minC, maxC) = circle.shape
    mask = np.zeros((M, N), dtype=int)
    minX = int(centX - minC / 2)
    maxX = int(centX + minC / 2)
    minY = int(centY - minC / 2)
    maxY = int(centY + maxC / 2)
    mask[minY:maxY, minX:maxX] = circle[:, :]

    FT = np.fft.fft2(field)
    FT = np.fft.fftshift(FT)
    filter = FT * mask
    #ft_holo = dis.intensity(filter, True)
    #ui.imageShow(ft_holo, 'FT hologram')

    min1X = int(N / 2 - (maxX - minX) / 2)
    max1X = int(N / 2 + (maxX - minX) / 2)
    min1Y = int(M / 2 - (maxY - minY) / 2)
    max1Y = int(M / 2 + (maxY - minY) / 2)

    crop_real = np.zeros((M, N))
    crop_real[min1Y:max1Y, min1X:max1X] = filter.real[minY:maxY, minX:maxX]

    crop_imag = np.zeros((M, N))
    crop_imag[min1Y:max1Y, min1X:max1X] = filter.imag[minY:maxY, minX:maxX]

    crop = crop_real + 1j * crop_imag
    ift = np.fft.ifftshift(crop)

    holoFilter = np.fft.ifft2(ift)
    field = holoFilter

    return field


def sfr(M, N, field, x1, x2, y1, y2):
    """
    # Function to create a spatial filter with a rectangle mask
    # Inputs:
    # M - Pixel size x
    # N - Pixel size y
    # field - The field to be filtered
    # x1 - Coordinate x1 for the rectangle (upper left corner)
    # y1 - Coordinate y1 for the rectangle (upper left corner)
    # x2 - Coordinate x2 for the rectangle (lower right corner)
    # y2 - Coordinate y2 for the rectangle (lower right corner)
    """
    mask = np.zeros((M, N))
    mask[y1:y2, x1:x2] = 1
    FT = np.fft.fft2(field)
    FT = np.fft.fftshift(FT)
    filter = FT * mask

    minX = int(N / 2 - (x2 - x1) / 2)
    maxX = int(N / 2 + (x2 - x1) / 2)
    minY = int(M / 2 - (y2 - y1) / 2)
    maxY = int(M / 2 + (y2 - y1) / 2)

    crop_real = np.zeros((M, N))
    crop_real[minY:maxY, minX:maxX] = filter.real[y1:y2, x1:x2]
    crop_imag = np.zeros((M, N))
    crop_imag[minY:maxY, minX:maxX] = filter.imag[y1:y2, x1:x2]

    crop = crop_real + 1j * crop_imag
    ift = np.fft.ifftshift(crop)

    holoFilter = np.fft.ifft2(ift)
    field = holoFilter

    return field