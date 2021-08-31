# -*- coding: utf-8 -*-
"""
Title-->            Example of Slightly Off-axis Synthetic reference (SOSR) quadrature phase shifting method.
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             01/05/2021
Last modified-->    31/08/2021
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Sample code to use the SOSR to recovered information of a Fresnel Lens from four pi/2-shifted DHM holograms.
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

import numpy as np

#import numericalPropagation.propagators as npr
import utilities
import utilities.display as ud
import phaseShifting.phaseShiftingMethods as phs


#pyDiffraction welcome message
utilities.salutation()

#Load four pi/2-shifted DHM holograms and automatically converts them into gray-scale images

Inp0 = utilities.imageRead('data/phaseShifting/I0_Fresnel_Lens.bmp')
Inp1 = utilities.imageRead('data/phaseShifting/I1_Fresnel_Lens.bmp')
Inp2 = utilities.imageRead('data/phaseShifting/I2_Fresnel_Lens.bmp')
Inp3 = utilities.imageRead('data/phaseShifting/I3_Fresnel_Lens.bmp')

#Background noise subtraction
Inp0 = Inp0-np.average(Inp0)
Inp1 = Inp1-np.average(Inp1)
Inp2 = Inp2-np.average(Inp2)
Inp3 = Inp3-np.average(Inp3)

#Parameters for the phase-shifting algorithm
#All units are calculated in meters
wavelength =  632.8e-9 #Illumination source wavelength
deltaX = deltaY = 6.9e-6 #Pixel pitch of camera sensor

#These variables must be carefully adjusted, see Documentation for further details
s = 1 #Size (in pixels) of the search grid
steps = 4 #Number of steps (along each direction) in the search grid

#Phase shifting computation via SOSR (Slighlty off-axis synthetic reference) method
comp_phase = phs.SOSR(Inp0, Inp1, Inp2, Inp3, True, wavelength, deltaX, deltaY, s, steps)

inten = ud.amplitude(comp_phase, False)
phase = ud.phase(comp_phase)

#Display an gray value image with the given title
utilities.imageShow(inten, 'Intensity of the sample')
utilities.imageShow(phase, 'Compensated phase of the sample')
