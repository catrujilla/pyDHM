# -*- coding: utf-8 -*-
"""
Title-->            Angular Spectrym Example
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             03/03/2019
Last modified-->    28/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Sample code to use the angular spectrum method to observe different Fresnel Zones
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

import numpy as np

import numericalPropagation
from utilities.display import show
import utilities

#pyDiffraction welcome message
utilities.salutation()

#Load an image and automatically converts it into a gray-scale image
hologram = utilities.imread('hol_testUSAF_Obj_8.33.bmp')

#Display an gray value image with the given title
imageShow (defocus_im, 'DHM hologram')

#Initially, the DHM hologram must be spatially-filtered
filt_holo = spatial_filtering (defocus_im-np.average(defocus_im))

#All units are calculated in meters
wavelength =  528e-9
deltaX = 4.687e-6

for i in range (-200,50,40):

    dist = i*1e-3
    #Propagating via AS to focus/defocus
    complexfield = angularSpectrum( filt_holo, dist, wavelength, deltaX, deltaX)

    #This function calculates the amplitude representation of a given complex field
    out = intensity(complexfield, False)

    #Display a gray-value image with the given title
    imageShow (out, 'Propagated image_'  + str(i) + ' mm' )