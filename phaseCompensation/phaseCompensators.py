"""
Title-->            Phase compensator package
Author-->           Carlos Trujillo, Ana Doblas and Raul Castaneda
Date-->             21/06/2021
Last modified-->    11/04/2022
Groups-->           University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
Abstract -->        Script for implementing four-phase compensators algorithms implemented in pyDHM:
                    Full ROI Search (FRS), Efficient ROI Search (ERS), the Cost Function ROI (CFS),
                    and the compensation No-tele (CNT).
Links-->            More information can be found on the website
                    https://github.com/catrujilla/pyDHM
information -->     This package uses the python library scipy for the CFS function. For the scipy
                    installation please visit the webpage https://scipy.org/install/
"""

# libraries
import numpy as np
import math
import phaseCompensation._init_ as pci
from math import pi
from scipy.optimize import minimize
import utilities.tools as tl

import utilities._init_ as ui
import utilities.display as dis


def FRS(inp, upper, wavelength, dx, dy, s=5, step=0.2):
    '''
    # Function to compensate phase maps of off-axis DHM via the full ROI search algorithm.
    # Inputs:
    # inp - The input intensity (captured) hologram
    # upper - Boolean variable defining whether the ROI will be that on the upper part of the DHM hologram spectrum.
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    # s = 2 and steps = 10
    '''

    # Retrieving the input shape
    inp = inp - np.average(inp)
    inp = np.array(inp)
    M, N = inp.shape

    # Creating a mesh_grid to operate in world-coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY

    # The spatial filtering process is executed
    print("Spatial filtering process started.....")
    holo_filter, indI, tmp = pci.spatialFiltering(inp, M, N, upper)
    print("Spatial filtering process finished.")

    # DC coordinates (center of DHM hologram Spectrum)
    fx_0, fy_0 = N / 2, M / 2

    # Coordintes of maxima in the DHM hologram spectrum
    fx_1, fy_1 = indI[1] + 1, indI[0] + 1  # +1 is the correction of the offset added by numpy?

    sum_max = 0
    arrayX = np.linspace(fx_1 - s, fx_1 + s, step)
    arrayY = np.linspace(fy_1 - s, fy_1 + s, step)

    k = (2 * pi) / wavelength

    print("Phase compensation started....")

    for fx_temp in arrayX:
        for fy_temp in arrayY:
            # print (fx_temp, fy_temp)
            theta_x = math.asin((fx_0 - fx_temp) * wavelength / (N * dx))
            theta_y = math.asin((fy_0 - fy_temp) * wavelength / (M * dy))
            ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
            # Compensation of the tilting angle for the off-axis acquisition
            reconstruction = holo_filter * ref_wave
            # module 2-pi phase retrieval
            phase = np.angle(reconstruction)
            # Thresholding process
            minVal = np.amin(phase)
            maxVal = np.amax(phase)
            phase_sca = (phase - minVal) / (maxVal - minVal)
            binary_phase = (phase_sca > 0.2)

            # plt.imshow(binary_phase, cmap='gray'), plt.title("binary_phase")  # image in gray scale
            # plt.show()  # show image

            # Applying the summation and thresholding metric
            sum = np.sum(np.sum(binary_phase))
            if sum > sum_max:
                x_max_out = fx_temp;
                y_max_out = fy_temp;
                sum_max = sum;
            # print (sum)

    # Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * dx))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * dy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    comp_phase = holo_filter * ref_wave

    print("Phase compensation finished.")
    # print(x_max_out, y_max_out)

    return comp_phase


def ERS(inp, upper, wavelength, dx, dy, s, step):
    '''
    # Function to compensate phase maps of off-axis DHM via the fast ROI search algorithm.
    # Inputs:
    # inp - The input intensity (captured) hologram
    # upper - Boolean variable defining whether the ROI will be that on the upper part of the DHM hologram spectrum.
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    # s = 5 and step = 0.2
    '''

    # Retrieving the input shape
    inp = inp - np.average(inp)
    inp = np.array(inp)
    M, N = inp.shape

    # Creating a mesh_grid to operate in world-coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY

    # The spatial filtering process is executed
    print("Spatial filtering process started.....")
    holo_filter, indI, tmp = pci.spatialFiltering(inp, M, N, upper)
    print("Spatial filtering process finished.")

    # out = phc.intensity(tmp, True)
    # plt.imshow(out, cmap='gray'), plt.title("holo_filterFreq")  # image in gray scale
    # plt.show()  # show image

    # Phase Compensation process

    # DC coordinates (center of DHM hologram Spectrum)
    fx_0, fy_0 = N / 2, M / 2

    # Coordintes of maxima in the DHM hologram spectrum
    fx_1, fy_1 = indI[1] + 1, indI[0] + 1  # +1 is the correction of the offset added by numpy?

    # print (fx_0, fy_0, fx_1, fy_1)

    print("Phase compensation started....")

    k = (2 * pi) / wavelength

    fin = 0
    sum_max = 0

    fx = fx_1
    fy = fy_1

    G_temp = s

    while fin == 0:

        sum_max = 0  # small number for the metric (thresholding)

        arrayY = np.linspace(int(10 * (fy - step * G_temp)), int(10 * (fy + step * G_temp)), int(10 * step))
        arrayX = np.linspace(int(10 * (fx - step * G_temp)), int(10 * (fx + step * G_temp)), int(10 * step))

        for fx_temp in arrayX:
            for fy_temp in arrayY:

                fx_tmp, fy_tmp = fx_temp / 10, fy_temp / 10

                # print (fx_tmp, fy_tmp)

                theta_x = math.asin((fx_0 - fx_tmp) * wavelength / (N * dx))
                theta_y = math.asin((fy_0 - fy_tmp) * wavelength / (M * dy))
                ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))

                # Compensation of the tilting angle for the off-axis acquisition
                reconstruction = holo_filter * ref_wave
                # module 2-pi phase retrieval
                phase = np.angle(reconstruction)
                # Thresholding process
                minVal = np.amin(phase)
                maxVal = np.amax(phase)
                phase_sca = (phase - minVal) / (maxVal - minVal)
                binary_phase = (phase_sca > 0.2)

                # plt.imshow(phase, cmap='gray'), plt.title("binary_phase")  # image in gray scale
                # plt.show()  # show image

                # Applying the summation and thresholding metric
                sum = np.sum(np.sum(binary_phase))
                if (sum > sum_max):
                    x_max_out = fx_tmp
                    y_max_out = fy_tmp
                    sum_max = sum
                # print (sum)

        G_temp = G_temp - 1

        if x_max_out == fx and y_max_out == fy:
            fin = 1;

        fx = x_max_out
        fy = y_max_out

    # print(x_max_out, y_max_out)

    # Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * dx))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * dy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    comp_phase = holo_filter * ref_wave

    print("Phase compensation finished.")

    return comp_phase


def CFS(inp, wavelength, dx, dy):
    '''
    # Function to compensate phase maps of image plane off-axis DHM via a cost-function search algorithm.
    # Inputs:
    # inp - The input intensity (captured) hologram
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    '''

    # Retrieving the input shape
    inp = np.array(inp)
    M, N = inp.shape

    # if M and N are not equals, cutting the hologram in a square shape
    if M > N:
        cut = round(int(M - N) / 2)
        inp = inp[cut:(M - cut):N, 0:N]
        M, N = inp.shape
        print('The hologram has been cut into a square shape of: [', M, 'x', N, ']')
    elif M < N:
        cut = round(int(N - M) / 2)
        inp = inp[0:M, cut:(N - cut)]
        M, N = inp.shape
        print('The hologram has been cut into a square shape of: [', M, 'x', N, ']')
    else:
        inp = inp
        M, N = inp.shape

    # Creating a mesh_grid to operate in world coordinates
    x = np.arange(0, M, 1)  # array x
    y = np.arange(0, N, 1)  # array y
    X, Y = np.meshgrid(x - (M / 2), y - (N / 2), indexing='xy')  # meshgrid XY
    fx_0 = M / 2
    fy_0 = N / 2
    k = (2 * math.pi) / wavelength

    # The spatial filtering process is executed
    print("Spatial filtering process started.....")
    holo_filter, fx_max, fy_max = pci.spatialFilteringCF(inp, M, N)
    print("Spatial filtering process finished.")

    # loading seeds
    seeds = [fx_max, fy_max]

    # minimization
    print("Minimization process started.....")
    step = 1.5
    res = minimize(pci.costFunction, seeds, args=(M, N, holo_filter, wavelength, dx, dy, X, Y, fx_0, fy_0, k),
                   method='TNC', bounds=((seeds[0] - step, seeds[0] + step), (seeds[1] - step, seeds[1] + step)), tol=1e-3)
    print("Minimization process finished.")
    fx_max = res.x[0]
    fy_max = res.x[1]
    print('x: ',fx_max)
    print('y: ',fy_max)
    fx_max = fx_max - step
    fy_max = fy_max - step

    # Best phase compensation
    print("Phase compensation started....")
    theta_x = math.asin((fx_0 - fx_max) * wavelength / (M * dx))
    theta_y = math.asin((fy_0 - fy_max) * wavelength / (N * dy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    comp_phase = holo_filter * ref_wave

    print("Phase compensation finished.")

    return comp_phase


def CNT(inp, wavelength, dx, dy, x1, x2, y1, y2, cur, s, step):
    '''
    # Function to compensate phase maps of image plane off-axis DHM, operating in non-telecentric regimen
    # Inputs:
    # inp - The input intensity (captured) hologram
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    '''

    # Retrieving the input shape
    inp = np.array(inp)
    M, N = inp.shape
    k = (2 * math.pi) / wavelength


    # Creating a mesh-grid to operate in world coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY


    # Fourier transform of the hologram
    FT = np.fft.fft2(inp)
    FT = np.fft.fftshift(FT)
    FT_display = dis.intensity(FT, True)


    # Normalization of the Fourier transform
    minVal = np.amin(FT_display)
    maxVal = np.amax(FT_display)
    FT_normalized = (FT_display - minVal) / (maxVal - minVal)
    binary_FT = (FT_normalized > 0.6)
    ui.imageShow(binary_FT, 'FT binirezed')


    # Filter to find the X_center and Y_center
    mask = np.zeros((M, N))
    mask[y1:y2, x1:x2] = 1
    filter = binary_FT * mask
    #ui.imageShow(filter, 'FT Filter')


    # X_center and Y_center
    Xmoment = 0
    Ymoment = 0
    Xaverage = 0
    Yaverage = 0
    for p in range(0, M):
        for q in range(0, N):
            Xmoment = Xmoment + filter[p, q] * q
            Xaverage = Xaverage + filter[p, q]
            Ymoment = Ymoment + filter[p, q] * p
            Yaverage = Yaverage + filter[p, q]
    Xcenter = Xmoment/Xaverage
    Ycenter = Ymoment/Yaverage


    # reference wave for the first compensation (global linear compensation)
    ThetaXM = math.asin((N / 2 - Xcenter) * wavelength / (M * dx))
    ThetaYM = math.asin((M / 2 - Ycenter) * wavelength / (N * dy))
    reference = np.exp(1j * k * (math.sin(ThetaXM) * X * dx + math.sin(ThetaYM) * Y * dy))


    # spatial filter
    mask = np.zeros((M, N))
    mask[y1 - 20:y2 + 20, x1 - 20:x2 + 20] = 1
    filter = FT * mask
    holo_filter = np.fft.ifftshift(filter)
    holo_filter = np.fft.ifft2(holo_filter)


    # First compensation
    comp_phase = holo_filter * reference


    # Binarization
    minVal = np.amin(comp_phase)
    maxVal = np.amax(comp_phase)
    phase_normalized = (comp_phase - minVal) / (maxVal - minVal)
    binary_phase = (phase_normalized > 0.6)
    ui.imageShow(binary_phase, 'Binarized phase')


    # X_center and Y_center
    Xmoment = 0
    Ymoment = 0
    Xaverage = 0
    Yaverage = 0
    for p in range(0, M):
        for q in range(0, N):
            Xmoment = Xmoment + binary_phase[p, q] * q
            Xaverage = Xaverage + binary_phase[p, q]
            Ymoment = Ymoment + binary_phase[p, q] * p
            Yaverage = Yaverage + binary_phase[p, q]
    Xcenter = Xmoment/Xaverage
    Ycenter = Ymoment/Yaverage

    p = input("Enter the pixel position on the x axis for the center of circular phase map ")
    q = input("Enter the pixel position on the y axis for the center of circular phase map  ")
    Xcenter = abs(M/2 - int(p))
    Ycenter = abs(N/2 - int(q))
    arrayXcenter = np.arange(Xcenter - s, Xcenter + s, step)
    arrayYcenter = np.arange(Ycenter - s, Ycenter + s, step)
    cont = 0
    sum_max = 0
    for Xcenter in arrayXcenter:
         for Ycenter in arrayYcenter:
            cont = cont + 1
            spheMAP = np.exp(-1j * (np.power((X - Xcenter), 2) + np.power((Y - Ycenter), 2)) / cur)
            phaseCompensate = comp_phase * spheMAP
            phaseCompensate = np.angle(phaseCompensate)
            #ui.imageShow(phaseCompensate, 'phaseCompensate')

            minVal = np.amin(phaseCompensate)
            maxVal = np.amax(phaseCompensate)
            phase_sca = (phaseCompensate - minVal) / (maxVal - minVal)
            binary_phase = (phase_sca > 0.2)
            #ui.imageShow(binary_phase, 'phaseCompensate')

            # Applying the summation and thresholding metric
            sum = np.sum(np.sum(binary_phase))
            if (sum > sum_max):
                xCenter_out = Xcenter
                yCenter_out = Ycenter
                sum_max = sum
                i = cont

    Xcenter = xCenter_out
    Ycenter = yCenter_out

    spheMAP = np.exp(-1j * (np.power((X - Xcenter), 2) + np.power((Y - Ycenter), 2)) / cur)
    phaseCompensate = comp_phase * spheMAP

    return phaseCompensate


