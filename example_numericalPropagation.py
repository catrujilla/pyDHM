"""
Title-->            Numerical propagation example.
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             12/20/2021
                    University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
Links-->          -

"""
import numpy as np
import utilities._init_ as ui
import utilities.display as dis
import numericalPropagation.propagators as pr
import utilities.tools as tl
import utilities.speckleMethods as speckle


# Angular spectrum

# Load the input plane
input = ui.imageRead('data/numericalPropagation samples/UofM-1-inv.jpg')
input = ui.imageRead('data/numericalPropagation samples/+3cm.tif')
ui.imageShow(input, 'input field')

# FT of the hologram
ft_holo = dis.FT(input)
ft_holo = dis.intensity(ft_holo, True)
ui.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = tl.sfc(input, 160, 303, 276)

# Numerical propagation using the angular spectrum
#output = pr.angularSpectrum(filter, 70000, 0.633, 6.9, 6.9)

# Display the output field
#intensity = dis.intensity(output, False)
#ui.imageShow(intensity, 'output field')

for z in range(-100, 80, 10):
    # Warning: Depending on the propagation distance step, this cycle can takke too much time.
    # Numerical propagation using the angular spectrum
    output = pr.angularSpectrum(filter, z, 0.000633, 0.0069, 0.0069)

    # Display the output field
    intensity = dis.intensity(output, False)
    ui.imageShow(intensity, 'output field')



# Fresnel transform and speckle reduction vu HM2F

# Load the input plane
input = ui.imageRead('data/numericalPropagation samples/horse.bmp')
ui.imageShow(input, 'input field')

# FT of the hologram
ft_holo = dis.FT(input)
ft_holo = dis.intensity(ft_holo, True)
ui.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = tl.sfr(input, 74, 500, 54, 450)

# Numerical propagation using the Fresnel transforms
output = pr.fresnel(filter, -450, 0.000633, 0.005000, 0.005000)

# Display the output field
amplitude = dis.amplitude(output, False)
ui.imageShow(amplitude, 'output field')

# HM2F to reduce the speckle
# denoise = speckle.HM2F(amplitude, 7, False, False)

#amplitude = dis.amplitude(denoise, False)
#ui.imageShow(amplitude, 'output denoise')



# Fresnel-Bluestein transform

# Load the input plane
input = ui.imageRead('data/numericalPropagation samples/die_1.jpg')
ui.imageShow(input, 'input field')

# FT of the hologram
ft_holo = dis.FT(input)
ft_holo = dis.intensity(ft_holo, True)
ui.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = tl.sfr(input, 280, 500, 150, 340)

# Numerical propagation using the Fresnel transforms
output = pr.bluestein(filter - np.average(filter), 30/100, 632.8e-9, 7e-6, 7e-6, 4.5e-5, 4.5e-5)

# Display the output field
amplitude = dis.amplitude(output, False)
ui.imageShow(amplitude, 'output field')


