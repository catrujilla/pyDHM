# -*- coding: utf-8 -*-
"""
Title-->            Fresnel Bluestein Example - Fourier method
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             09/07/2021
Last modified-->    28/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Sample code to use the fresnel Bluestein method to numerically reconstruct an off-axis hologram
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

import numpy as np

import numericalPropagation
from utilities.display import imageShow
import utilities

#pyDiffraction welcome message
utilities.salutation()

#Load an image and automatically converts it into a gray-scale image
hologram = utilities.imageRead('data/HoloMedal_3.tif')

#Display an gray value image with the given title
imageShow (hologram, 'Hologram')

#This function calculates the propagated complex field of an input transmitance 
#(in this case this function is used to numerically reconstruct an hologram)
i = 380
#complexfield = numericalPropagation.fresnel( (hologram - np.average(hologram)), 632.8e-9, i/100, 11e-6, 11e-6)

for i in range (250,350,10):
    
    complexfield = numericalPropagation.bluestein( (hologram - np.average(hologram)), i/100, 532, 4.5e-6, 4.5e-6, 4.5e-6, 4.5e-6)

    #This function calculates the amplitude representatin of a given complex field
    #(The second parameter allows to determine if a log operation is performed)
    out = utilities.display.amplitude(complexfield, True)

    #Display an gray value image with the given title
    imageShow (out, 'Amplitude recontruction - log display ' + str(i) + ' cms' )

#Save this data into a image file in the disk
#utilities.imageSave ('Amplitude recontruction - log display ' + str(i) + ' cms.bmp' , out)
utilities.imageSave ('recontruction.bmp' , out)
