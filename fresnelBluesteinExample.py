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

import numericalPropagation.propagators as npr
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
deltaX_out = deltaY_out = 8.5e-5 #output pixel pitch
 
#This function calculates the propagated complex field of an input transmitance 
#(in this case this function is used to numerically reconstruct an hologram with scaled output coordinates)
complexfield = npr.bluestein( (hologram - np.average(hologram)), z/100, wavelength, deltaX_in, deltaY_in, deltaX_out, deltaY_out)
out = ud.amplitude(complexfield, False)

#Display an gray value image with the given title
utilities.imageShow(out, 'Amplitude recontruction - log display ' + str(z) + ' cms' )

