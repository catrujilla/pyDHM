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
from scipy import ndimage
from matplotlib import pyplot as plt


def HM2F(inp, kernel, figures, plots):
    if kernel % 2 == 0:
        print('Kernel size must be a odd number')
        exit()

    mean_image = inp
    cont = 1
    for i in range(3, kernel + 2, 2):
        filter = ndimage.median_filter(inp, i, mode='constant', cval=0)
        mean_image = (mean_image + filter) / 2
        imageDenoise = mean_image

    return imageDenoise


