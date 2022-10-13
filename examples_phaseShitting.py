"""
Title-->            Phase-shifting reconstruction examples.
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             11/04/2022
Groups-->           University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
e-mail-->           catrujilla@eafit.edu.co and adoblas@memphis.edu
Links-->            More information about the library can be found on the website
                    https://github.com/catrujilla/pyDHM         -
"""

from pyDHM import utilities
from pyDHM import phaseShifting

# Phase shifting using the SOSR approach
"""
print ("SOSR example")
# Load the holograms
inp0 = utilities.imageRead('data/phase-shifting samples/cut_samples/I0_Fresnel_Lens.bmp')
inp1 = utilities.imageRead('data/phase-shifting samples/cut_samples/I1_Fresnel_Lens.bmp')
inp2 = utilities.imageRead('data/phase-shifting samples/cut_samples/I2_Fresnel_Lens.bmp')
inp3 = utilities.imageRead('data/phase-shifting samples/cut_samples/I3_Fresnel_Lens.bmp')

# Phase shifting using the SOSR approach
output = phaseShifting.SOSR(inp0, inp1, inp2, inp3, True, 632.8e-9, 6.9e-6, 6.9e-6, 1, 4)

# Display the amplitude reconstruction
amplitude = utilities.amplitude(output, False)
utilities.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
"""

# Phase shifting via ps5-inline

print ("ps5-inline example")

inp0 = utilities.imageRead('data/phase-shifting samples/PS5/h1_PS5.png')
inp1 = utilities.imageRead('data/phase-shifting samples/PS5/h2_PS5.png')
inp2 = utilities.imageRead('data/phase-shifting samples/PS5/h3_PS5.png')
inp3 = utilities.imageRead('data/phase-shifting samples/PS5/h4_PS5.png')
inp4 = utilities.imageRead('data/phase-shifting samples/PS5/h5_PS5.png')

# Phase shifting via ps5 ()
output = phaseShifting.PS5(inp0, inp1, inp2, inp3, inp4)
phase = utilities.phase(output)

# Display the phase reconstruction
utilities.imageShow(phase, 'Phase reconstruction')


# Phase shifting via ps4-inline

print ("ps4-inline example")
inp0 = utilities.imageRead('data/phase-shifting samples/PS4/h1-PS4.png')
inp1 = utilities.imageRead('data/phase-shifting samples/PS4/h2-PS4.png')
inp2 = utilities.imageRead('data/phase-shifting samples/PS4/h3-PS4.png')
inp3 = utilities.imageRead('data/phase-shifting samples/PS4/h4-PS4.png')

# Phase shifting via ps4 ()
output = phaseShifting.PS4(inp0, inp1, inp2, inp3)
phase = utilities.phase(output)

# Display the phase reconstruction
utilities.imageShow(phase, 'Phase reconstruction')


# Phase shifting via ps3-inline

print ("ps3-inline example")
inp0 = utilities.imageRead('data/phase-shifting samples/PS3/h1_PS3.png')
inp1 = utilities.imageRead('data/phase-shifting samples/PS3/h2_PS3.png')
inp2 = utilities.imageRead('data/phase-shifting samples/PS3/h3_PS3.png')

# Phase shifting via ps3 ()
output = phaseShifting.PS3(inp0, inp1, inp2)
phase = utilities.phase(output)

# Display the phase reconstruction
utilities.imageShow(phase, 'Phase reconstruction')


# three blind raw Holograms
"""
print ("BPS3 example")
# Load the holograms
inp0 = utilities.imageRead('data/phase-shifting samples/usaf_1_cut.jpg')
inp1 = utilities.imageRead('data/phase-shifting samples/usaf_2_cut.jpg')
inp2 = utilities.imageRead('data/phase-shifting samples/usaf_3_cut.jpg')

utilities.imageShow(inp0, 'Hologram 1')
utilities.imageShow(inp1, 'Hologram 2')
utilities.imageShow(inp2, 'Hologram 3')

# Phase shifting via BPS3 (three holograms slightly off-axis)
output = phaseShifting.BPS3(inp2, inp1, inp0, 0.532, 2.4, 2.4)

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')


# two blind raw Holograms

print ("BPS2 example")
# Load the holograms
#inp0 = utilities.imageRead('data/phase-shifting samples/Neuron_1.jpg')
#inp1 = utilities.imageRead('data/phase-shifting samples/Neuron_3.jpg')
inp0 = utilities.imageRead('data/phase-shifting samples/usaf_1_cut.jpg')
inp1 = utilities.imageRead('data/phase-shifting samples/usaf_2_cut.jpg')

# Phase shifting via BPS2 (two holograms slighlty off-axis)
output = phaseShifting.BPS2(inp1, inp0, 0.532, 2.4, 2.4)

# Display the phase reconstruction
phase = utilities.phase(output)
utilities.imageShow(phase, 'Phase reconstruction')
"""