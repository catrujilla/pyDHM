"""
Title-->            Slightly off-axis numerical phase compensation examples.
Author-->           Carlos Trujillo, Ana Doblas, and Raul Castaneda
Date-->             11/04/2022
Groups-->           University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
e-mail-->           catrujilla@eafit.edu.co and adoblas@memphis.edu
Links-->            More information about the library can be found on the website
                    https://github.com/catrujilla/pyDHM
"""
import numpy as np
from pyDHM import utilities
from pyDHM import phaseCompensation
from pyDHM import numericalPropagation

from timeit import default_timer as timer

# FRS example using library
'''
print ("FRS example")
# Load the hologram
hologram = utilities.imageRead('data/compensation samples/fly 20x.bmp')
utilities.imageShow(hologram, 'Hologram')

start = timer()	#Start to count time

# Numerical compensation using the FRS approach
output = phaseCompensation.FRS(hologram, True, 0.532, 2.4, 2.4, 2, 10)

print("With CPU-only processing:", timer()-start) #Time for fft_shift execution

# Display the amplitude reconstruction
amplitude = utilities.amplitude(output, False)
utilities.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
'''

# ERS example
'''
print ("ERS example")
# Load the hologram
#hologram = utilities.imageRead('data/compensation samples/star-target.png')
hologram = utilities.imageRead('data/compensation samples/fly 20x.bmp')
utilities.imageShow(hologram, 'Hologram')

# Numerical compensation using the ERS approach
#output = phaseCompensation.ERS(hologram, True, 0.532, 2.6, 2.6, 5, 0.2)
output = phaseCompensation.ERS(hologram, True, 0.532, 2.4, 2.4, 5, 0.2)

# Display the amplitude reconstruction
amplitude = utilities.amplitude(output, False)
utilities.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
'''

# CFS example
'''
print ("CFS example")
# Load the hologram
#hologram = utilities.imageRead('data/compensation samples/hologram EAFIT.jpg')
hologram = utilities.imageRead('data/compensation samples/fly 20x.bmp')
utilities.imageShow(hologram, 'Hologram')

# Numerical compensation using the CFS approach
output = phaseCompensation.CFS(hologram, 0.532, 2.4, 2.4)

# Display the amplitude reconstruction
amplitude = utilities.amplitude(output, False)
utilities.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
'''

# CNT example

print ("CNT example")
# Load the hologram
hologram = utilities.imageRead('data/compensation samples/holonotele_fly-resize.png')
# utilities.imageShow(hologram, 'Hologram')

# FT of the hologram
ft_holo = utilities.FT(hologram)
ft_holo = utilities.intensity(ft_holo, True)
#utilities.imageShow(ft_holo, 'FT hologram')

# Numerical compensation using the CNT approach
output = phaseCompensation.CNT(hologram, 0.633, 7, 7, spatialFilter='sfmr')

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')


# for the propagation
'''
for i in range(0, 10000, 1000):
    dist = i
    # Propagating via angular spectrum to focus the information of the sample
    complexfield = numericalPropagation.angularSpectrum(output, dist, 0.532, 5.4, 5.4)

    # This function calculates the amplitude representation of a given complex field
    out = utilities.amplitude(complexfield, False)

    # Display a gray-value image with the given title
    utilities.imageShow(out, 'Propagated amplitude image_' + str(i) + ' mm')

    # This function calculates the amplitude representation of a given complex field
    out = utilities.phase(complexfield)
    utilities.imageShow(out, 'Propagated phase image_' + str(i) + ' mm')
'''
