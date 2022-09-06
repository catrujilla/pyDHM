
#Fig. 3
# import packages
from pyDHM import utilities
from pyDHM import phaseShifting

# Load the holograms
inp0 = utilities.imageRead('holo1.png')
inp1 = utilities.imageRead('holo2.png')
inp2 = utilities.imageRead('holo3.png')
inp3 = utilities.imageRead('holo4.png')
inp4 = utilities.imageRead('holo5.png')

# Phase shifting using the PS strategies
output = phaseShifting.PS5(inp0, inp1, inp2, inp3, inp4)
output = phaseShifting.PS4(inp0, inp1, inp2, inp3)
output = phaseShifting.PS3(inp0, inp1, inp2)

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(output, 'Phase reconstruction') 


#Fig. 4
# import packages 
from pyDHM import utilities
from pyDHM import phaseShifting 

# Load the holograms
inp0 = utilities.imageRead('holo1.png')
inp1 = utilities.imageRead('holo2.png')
inp2 = utilities.imageRead('holo3.png')
inp3 = utilities.imageRead('holo4.png')

# Phase shifting via SOSR, BPS3 or BPS2
output = phaseShifting.SOSR(inp0, inp1, inp2, inp3, 633e-9, 6.9e-6, 6.9e-6, 1, 4)
output = phaseShifting.BPS3(inp0, inp1, inp2, 0.532, 2.9, 2.9)
output = phaseShifting.BPS2(inp0, inp1, 0.532, 2.9, 2.9)

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')


'''
#Fig. 5
# import packages 
from pyDHM import utilities
from pyDHM import phaseCompensation

# Load the hologram
inp = utilities.imageRead(‘hologram.jpg’)

# Compensation using the available approaches
output = phaseCompensation.FRS(inp, True, 0.633, 6.9, 6.9, 2, 10)
output = phaseCompensation.ERS(inp, True, 0.633, 6.9, 6.9, 5, 0.2)
output = phaseCompensation.CFS(inp, 0.532, 2.6, 2.6)

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
'''

'''
#Fig. 6 
# import packages 
from pyDHM import utilities
from pyDHM import phaseCompensation

# Load the hologram
inp = utilities.imageRead('hologram.tif')

# FT of the hologram
ft_holo = utilities.FT(inp)
ft_holo = utilities.intensity(ft_holo, True)
utilities.imageShow(ft_holo, 'FT hologram')

# Numerical compensation using the CNT approach
output = phaseCompensation.CNT(inp, 0.633, 6.9, 6.9, spatialFilter='sfmr')

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
'''

'''
#Fig 7
# import packages 
from pyDHM import utilities
from pyDHM import numericalPropagation

# Load the hologram
hologram =  utilities.imageRead('hologram.jpg')
utilities.imageShow(hologram, 'Hologram’)

# FT of the hologram
ft_holo = utilities.FT(hologram)
ft_holo = utilities.intensity(ft_holo, True)
utilities.imageShow(ft_holo, 'FT hologram’)

# Circular spatial filter
filter = utilities.sfc(hologram, 160, 303, 276)

# Numerical propagation using the angular spectrum
output = numericalPropagation.angularSpectrum(filter, z, 0.633, 6.9, 6.9)

# Display the output field 
intensity = utilities.intensity(output, False)
utilities.imageShow(intensity,‘ouput field’)
'''

'''
#Fig. 8
# import packages 
from pyDHM import utilities
from pyDHM import numericalPropagation

# Load the hologram
input = utilities.imageRead(hologram.bmp')
utilities.imageShow(input, 'input field’)

# FT of the hologram
ft_holo = dis.FT(input)
ft_holo = dis.intensity(ft_holo, True)
utilities.imageShow(ft_holo, 'FT hologram')

# Spatial filter
filter = utilities.sfr(input, 280, 500, 150, 340)

# Numerical propagation using the Fresnel transforms
output = numericalPropagation.fresnel(input, -450, 0.000633, 0.005, 0.005)
output =  numericalPropagation.bluestein(input, 0.03, 633e-9, 7.4e-6, 7.4e-6, 1.48e-5, 1.48e-5)                                                                                                    1.48e-5)

# Display the output field
amplitude = utilities.amplitude(output, True)
utilities.imageShow(amplitude, output fieldn')
'''
