# -*- coding: utf-8 -*-
"""
Title-->            Fresnel Example - Fourier method
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             03/03/2019
Last modified-->    16/07/2020
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Sample code to use the fresnel method to numerically reconstruct an off-axis hologram
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

import cv2
import numpy as np
from matplotlib import pyplot as plt
from numericalPropagation.fresnel import fresnel
from numericalPropagation.display import amplitude

namefile = 'HoloMedal_3.tif'
hologram = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
plt.imshow(hologram, cmap='gray'), plt.title('hologram')  # image in gray scale
plt.show()  # show hologram

print (np.average(hologram))

#This function calculates the propagated complex field of an input transmitance 
#(in this case this function is used to numerically reconstruct an hologram)
complexfield = fresnel( (hologram - np.average(hologram)), 532e-9, 3.480, 4.5e-6, 4.5e-6)

#This function calculates the amplitude representatin of a given complex field
#(The second parameter allows to determine if a log operation is performed)
out =  amplitude(complexfield, True)

plt.imshow(out, cmap='gray'), plt.title('Amplitude recontruction - log display')  # image in gray scale
plt.show()  # show hologram
