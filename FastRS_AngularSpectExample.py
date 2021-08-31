# -*- coding: utf-8 -*-
"""
Title-->            Fast ROI search with Angular Spectrum Example
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             17/06/2021
Last modified-->    31/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Sample code to use the angular spectrum method to propagate a complex wavefield compensated in phase with the FastRS algorithm.
Links-->          - https://github.com/catrujilla/pyDHM
"""

import numpy as np

import numericalPropagation.propagators as npr
import utilities
import utilities.display as ud
import phaseCompensation.phaseCompensators as phc


#pyDiffraction welcome message
utilities.salutation()

#Load an image and automatically converts it into a gray-scale image
defocus_hol = utilities.imageRead('data/Nymphadea2-1.jpg')

#Parameters for the retrieval of the compensated phase sample information
#All units are calculated in meters
wavelength =  528e-9 #Illumination source wavelength
deltaX = deltaY = 4.687e-6 #Pixel pitch of camera sensor

#These variables must be carefully adjusted, see Documentation for further details
s = 5 #The initial number of steps in each direction of the search
step = 0.2 #The step in the grid search

#Phase compensation computation via FRS (Full ROI search)
comp_phase = phc.FastRS(defocus_hol-np.average(defocus_hol), True, wavelength, deltaX, deltaY, s, step)

inten = ud.amplitude(comp_phase, False)
phase = ud.phase(comp_phase)

#Display an gray value image with the given title
utilities.imageShow(inten, 'Intensity of the sample')
utilities.imageShow(phase, 'Compensated phase of the sample')

'''
#If the sample is out-of-focus (no image-plane acqusition), we can focus by propagating.
#Propagating via angular spectrum method implemented in pyDHM.

for i in range (-100,100,20):

    dist = i*1e-3
    #Propagating via angular spectrum to focus the information of the sample
    complexfield = npr.angularSpectrum( comp_phase, dist, wavelength, deltaX, deltaX)

    #This function calculates the amplitude representation of a given complex field
    out = ud.amplitude(complexfield, False)
    
    #Display a gray-value image with the given title
    utilities.imageShow(out, 'Propagated image_'  + str(i) + ' mm' )
    
    #This function calculates the amplitude representation of a given complex field
    out = ud.phase(complexfield)
    utilities.imageShow(out, 'Propagated phase image_'  + str(i) + ' mm' )

'''