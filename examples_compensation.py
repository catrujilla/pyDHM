"""
Title-->            Off-axis numerical reconstruction example.
Author-->           Carlos Trujillo, Ana Doblas, and Raul Castaneda
Date-->             11/04/2022
Groups-->           University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
e-mail-->           catrujilla@eafit.edu.co and adoblas@memphis.edu
Links-->            More information about the library can be found on the website
                    https://github.com/catrujilla/pyDHM
"""
import numpy as np

import utilities._init_ as ui
import utilities.display as dis
import utilities.tools as tl
import phaseCompensation.phaseCompensators as pc
import numericalPropagation.propagators as pr


# FRS example using library
'''
print ("FRS example")
# Load the hologram
hologram = ui.imageRead('data/compensation samples/fly 20x.bmp')
ui.imageShow(hologram, 'Hologram')

# Numerical compensation using the FRS approach
output = pc.FRS(hologram, True, 0.532, 2.4, 2.4, 2, 10)

# Display the amplitude reconstruction
amplitude = dis.amplitude(output, False)
ui.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')
'''

# ERS example
'''
print ("ERS example")
# Load the hologram
#hologram = ui.imageRead('data/compensation samples/star-target.png')
hologram = ui.imageRead('data/compensation samples/fly 20x.bmp')
ui.imageShow(hologram, 'Hologram')

# Numerical compensation using the ERS approach
#output = pc.ERS(hologram, True, 0.532, 2.6, 2.6, 5, 0.2)
output = pc.ERS(hologram, True, 0.532, 2.4, 2.4, 5, 0.2)

# Display the amplitude reconstruction
amplitude = dis.amplitude(output, False)
ui.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')
'''

# CFS example
'''
print ("CFS example")
# Load the hologram
#hologram = ui.imageRead('data/compensation samples/hologram EAFIT.jpg')
hologram = ui.imageRead('data/compensation samples/fly 20x.bmp')
ui.imageShow(hologram, 'Hologram')

# Numerical compensation using the CFS approach
output = pc.CFS(hologram, 0.532, 2.4, 2.4)

# Display the amplitude reconstruction
amplitude = dis.amplitude(output, False)
ui.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')
'''

'''
# CNT example
print ("CNT example")
# Load the hologram
hologram = ui.imageRead('data/compensation samples/holonotele_fly-resize.png')
# ui.imageShow(hologram, 'Hologram')

# FT of the hologram
ft_holo = dis.FT(hologram)
ft_holo = dis.intensity(ft_holo, True)
#ui.imageShow(ft_holo, 'FT hologram')


# Numerical compensation using the CNT approach
output = pc.CNT(hologram, 0.633, 7, 7, 200, 287, 180, 267, 4200, 5, 1)


# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')

#phase_unw = np.apply_over_axes(np.unwrap, phase, np.arange(len(phase.shape)))
#phase_unw = np.unwrap(phase, discont=None, axis=- 1, period=3.143185307179586)
#ui.imageShow(phase_unw, 'phase_unw')
'''

# for the propagation
'''
for i in range(0, 10000, 1000):
    dist = i
    # Propagating via angular spectrum to focus the information of the sample
    complexfield = pr.angularSpectrum(output, dist, 0.532, 5.4, 5.4)

    # This function calculates the amplitude representation of a given complex field
    out = dis.amplitude(complexfield, False)

    # Display a gray-value image with the given title
    ui.imageShow(out, 'Propagated amplitude image_' + str(i) + ' mm')

    # This function calculates the amplitude representation of a given complex field
    out = dis.phase(complexfield)
    ui.imageShow(out, 'Propagated phase image_' + str(i) + ' mm')
'''
