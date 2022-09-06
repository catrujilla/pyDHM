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
output = phaseCompensation.CNT(inp, 0.633, 6.9, 6.9, 200, 287, 180, 267)

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
