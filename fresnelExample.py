# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 09:34:31 2020

@author: Calelo
"""

import cv2
import numpy as np
from matplotlib import pyplot as plt
from numericalPropagation.fresnel import fresnel
from numericalPropagation.fresnel import fr

namefile = 'HoloMedal_3.tif'
hologram = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
plt.imshow(hologram, cmap='gray'), plt.title('hologram')  # image in gray scale
plt.show()  # show hologram

print (np.average(hologram))

complexAmplitude = fresnel( (hologram - np.average(hologram)), 532e-9, 3.480, 4.5e-6, 4.5e-6)
amplitude = (np.abs(complexAmplitude))
amplitude = 20 * np.log(amplitude)  # logaritm scale FFT

plt.imshow(amplitude, cmap='gray'), plt.title('Amplitude recontruction - log display')  # image in gray scale
plt.show()  # show hologram

#def read_sample():
#    #namefile = 'mask.jpg'
#    namefile = 'die_1.jpg'
#    hologram = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
#    plt.imshow(hologram, cmap='gray'), plt.title('hologram')  # image in gray scale
#    plt.show()  # show hologram
#    return hologram
#
#def main():
#    hologram = read_sample()
#    # mask
#    intensitie = fr(hologram, 632.8, 1.05, 5.2, 5.2)
#
#    plt.imshow(intensitie, cmap='gray'), plt.title('phase recontruction')  # image in gray scale
#    plt.show()  # show hologram
#
#main()