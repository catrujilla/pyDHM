"""
Title-->            off-axis numerical reconstruction example.
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             12/20/2021
                    University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
Links-->          -

"""

import utilities._init_ as ui
import phaseShifting.phaseShifting as ps
import utilities.display as dis


# Phase shifting using the SOSR approach
'''
# Load the holograms
inp0 = ui.imageRead('data/phase-shifting samples/cut_samples/I0_Fresnel_Lens.bmp')
inp1 = ui.imageRead('data/phase-shifting samples/cut_samples/I1_Fresnel_Lens.bmp')
inp2 = ui.imageRead('data/phase-shifting samples/cut_samples/I2_Fresnel_Lens.bmp')
inp3 = ui.imageRead('data/phase-shifting samples/cut_samples/I3_Fresnel_Lens.bmp')

# Phase shifting using the SOSR approach
output = ps.SOSR(inp0, inp1, inp2, inp3, True, 632.8e-9, 6.9e-6, 6.9e-6, 1, 4)

# Display the amplitude reconstruction
amplitude = dis.amplitude(output, False)
ui.imageShow(amplitude, 'Amplitude reconstruction')

# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')
#'''


# Phase shifting via ps5-inline ()
'''
inp0 = ui.imageRead('data/phase-shifting samples/PS5/h1_PS5.png')
inp1 = ui.imageRead('data/phase-shifting samples/PS5/h2_PS5.png')
inp2 = ui.imageRead('data/phase-shifting samples/PS5/h3_PS5.png')
inp3 = ui.imageRead('data/phase-shifting samples/PS5/h4_PS5.png')
inp4 = ui.imageRead('data/phase-shifting samples/PS5/h5_PS5.png')

# Phase shifting via ps5 ()
output = ps.PS5(inp0, inp1, inp2, inp3, inp4)

# Display the phase reconstruction
ui.imageShow(output, 'Phase reconstruction')
'''


# Phase shifting via ps4-inline ()
'''
inp0 = ui.imageRead('data/phase-shifting samples/PS4/h1-PS4.png')
inp1 = ui.imageRead('data/phase-shifting samples/PS4/h2-PS4.png')
inp2 = ui.imageRead('data/phase-shifting samples/PS4/h3-PS4.png')
inp3 = ui.imageRead('data/phase-shifting samples/PS4/h4-PS4.png')

# Phase shifting via ps4 ()
output = ps.PS4(inp0, inp1, inp2, inp3)

# Display the phase reconstruction
ui.imageShow(output, 'Phase reconstruction')
'''


# Phase shifting via ps3-inline ()
'''
inp0 = ui.imageRead('data/phase-shifting samples/PS3/h1_PS3.png')
inp1 = ui.imageRead('data/phase-shifting samples/PS3/h2_PS3.png')
inp2 = ui.imageRead('data/phase-shifting samples/PS3/h3_PS3.png')

# Phase shifting via ps3 ()
output = ps.PS3(inp0, inp1, inp2)

# Display the phase reconstruction
ui.imageShow(output, 'Phase reconstruction')
'''


# three blind raw Holograms
'''
# Load the holograms
inp0 = ui.imageRead('data/phase-shifting samples/usaf_1_cut.jpg')
inp1 = ui.imageRead('data/phase-shifting samples/usaf_2_cut.jpg')
inp2 = ui.imageRead('data/phase-shifting samples/usaf_3_cut.jpg')

#ui.imageShow(inp0, 'Hologram 1')
#ui.imageShow(inp0, 'Hologram 2')
#ui.imageShow(inp0, 'Hologram 3')

# Phase shifting via BPS3 (three holograms slightly off-axis)
output = ps.BPS3(inp2, inp1, inp0, 0.532, 2.4, 2.4)

# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')
'''


# two blind raw Holograms
'''
# Load the holograms
inp0 = ui.imageRead('data/phase-shifting samples/Neuron_1.tif')
inp1 = ui.imageRead('data/phase-shifting samples/Neuron_3.tif')
inp0 = ui.imageRead('data/phase-shifting samples/usaf_1_cut.jpg')
inp1 = ui.imageRead('data/phase-shifting samples/usaf_2_cut.jpg')

# Phase shifting via BPS2 (two holograms slighlty off-axis)
output = ps.BPS2(inp1, inp0, 0.532, 2.4, 2.4)

# Display the phase reconstruction
phase = dis.phase(output)
ui.imageShow(phase, 'Phase reconstruction')
'''