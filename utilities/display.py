# -*- coding: utf-8 -*-
"""
Title-->            Display options utility script
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             03/03/2019
Last modified-->    16/07/2020
                    University of Memphis
                    Optical Imaging Research lab (OIRL)
                    EAFIT University
                    Applied Optics Group
Abstract -->        Script that implements the different methods to render the resulting complex field data
Links-->          - https://unal-optodigital.github.io/JDiffraction/
"""

import numpy as np

# Function to calcule the amplitude representation of a given complex field
# Inputs:
# inp - The input complex field
# log - boolean variable to determine if a log representation is applied

def amplitude (inp, log):
    
    out = np.abs(inp)
    
    if log == True:
        out = 20 * log(out)
	
    return out

