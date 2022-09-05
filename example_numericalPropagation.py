"""
Title-->            Numerical propagation examples.
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             05/09/2022
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Script that implements the different methods to propagate the resulting complex field data
Links-->          - https://github.com/catrujilla/pyDHM

"""
from pyDHM import utilities
from pyDHM import numericalPropagation
import numpy as np

# Angular spectrum
print ("Angular spectrum example")
# Load the input plane
inp = utilities.imageRead('data/numericalPropagation samples/UofM-1-inv.jpg')
inp = utilities.imageRead('data/numericalPropagation samples/+3cm.tif')
utilities.imageShow(inp, 'input field')

# FT of the hologram
ft_holo = utilities.FT(inp)
ft_holo = utilities.intensity(ft_holo, True)
utilities.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = utilities.sfmr(inp, True)

# Numerical propagation using the angular spectrum
output = numericalPropagation.angularSpectrum(filter, 70000, 0.633, 6.9, 6.9)

# Display the output field
intensity = utilities.intensity(output, False)
utilities.imageShow(intensity, 'focused output field')

for z in range(-100, 80, 10):
    # Warning: Depending on the propagation distance step, this cycle can takke too much time.
    # Numerical propagation using the angular spectrum
    output = numericalPropagation.angularSpectrum(filter, z, 0.000633, 0.0069, 0.0069)

    # Display the output field
    intensity = utilities.intensity(output, False)
    utilities.imageShow(intensity, 'output field ' + str(z))


print ("Fresnel transform example")
# Fresnel transform and speckle reduction vu HM2F

# Load the input plane
inp = utilities.imageRead('data/numericalPropagation samples/horse.bmp')
utilities.imageShow(inp, 'input field')

# FT of the hologram
ft_holo = utilities.FT(inp)
ft_holo = utilities.intensity(ft_holo, True)
utilities.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = utilities.sfmr(inp, True)

# Numerical propagation using the Fresnel transforms
output = numericalPropagation.fresnel(filter, -450, 0.000633, 0.005000, 0.005000)

# Display the output field
amplitude = utilities.amplitude(output, False)
utilities.imageShow(amplitude, 'output field')

# HM2F to reduce the speckle
denoise = utilities.HM2F(amplitude, 7, False, False)

#amplitude = utilities.amplitude(denoise, False)
utilities.imageShow(amplitude, 'output field denoised')


# Fresnel-Bluestein transform
print ("Fresnel-Bluestein transform example")
# Load the input plane
inp = utilities.imageRead('data/numericalPropagation samples/die_1.jpg')
utilities.imageShow(inp, 'input field')

# FT of the hologram
ft_holo = utilities.FT(inp)
ft_holo = utilities.intensity(ft_holo, True)
utilities.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = utilities.sfmr(inp, True)

# Numerical propagation using the Fresnel transforms
output = numericalPropagation.bluestein(filter - np.average(filter), 30/100, 632.8e-9, 7e-6, 7e-6, 4.5e-5, 4.5e-5)

# Display the output field
amplitude = utilities.amplitude(output, False)
utilities.imageShow(amplitude, 'output field')

