# libraries
import numpy as np
from math import pi

def info():
    
    print ("The available phase shifting algorithms are slightly off-axis synthetic reference (SOSR), blind shift (BS) and the traditinal full in-line 3-step, 4-step and 5-step algorithms")
    
    return
    
# Function to calcule the intensity representation of a given complex field
# Inputs:
# inp - The input complex field
# log - boolean variable to determine if a log representation is applied
def intensity (inp, log):
    
    out = np.abs(inp)
    out = out*out
    
    if log == True:
        out = 20 * np.log(out)
	
    return out