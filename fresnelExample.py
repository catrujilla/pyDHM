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

import numpy as np

import numericalPropagation as npr
#from utilities.display import imageShow
import utilities
import utilities.display as ud

#pyDiffraction welcome message
utilities.salutation()

#Load an image and automatically converts it into a gray-scale image
hologram = utilities.imageRead('data/die_1.jpg')

#Display an gray value image with the given title
utilities.imageShow(hologram, 'Hologram')

z = 78 #reconstruction (propagation) distance in cm
wavelength =  632.8e-9 #wavelength of illumination
deltaX_in = deltaY_in = 11e-6 #input pixel pitch

#This function calculates the propagated complex field of an input transmitance 
#(in this case this function is used to numerically reconstruct an hologram)
complexfield = npr.fresnel( (hologram - np.average(hologram)), z/100, wavelength, deltaX_in, deltaY_in)

#This function calculates the amplitude representatin of a given complex field
#(The second parameter allows to determine if a log operation is performed)
out = ud.amplitude(complexfield, False)

#Display an gray value image with the given title
utilities.imageShow(out, 'Amplitude recontruction - log display ' + str(z) + ' cms' )


