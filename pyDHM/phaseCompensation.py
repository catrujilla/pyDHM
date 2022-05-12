"""
Title-->            Phase compensators package
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
from matplotlib import pyplot as plt
from math import pi
from scipy.optimize import minimize

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
    holo_filter, indI, tmp = spatialFiltering(inp, M, N, upper)
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
    holo_filter, indI, tmp = spatialFiltering(inp, M, N, upper)
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
    holo_filter, fx_max, fy_max = spatialFilteringCF(inp, M, N)
    print("Spatial filtering process finished.")

    # loading seeds
    seeds = [fx_max, fy_max]

    # minimization
    print("Minimization process started.....")
    step = 1.5
    res = minimize(costFunction, seeds, args=(M, N, holo_filter, wavelength, dx, dy, X, Y, fx_0, fy_0, k),
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
    FT_display = intensity(FT, True)


    # Normalization of the Fourier transform
    minVal = np.amin(FT_display)
    maxVal = np.amax(FT_display)
    FT_normalized = (FT_display - minVal) / (maxVal - minVal)
    binary_FT = (FT_normalized > 0.6)
    imageShow(binary_FT, 'FT binirezed')


    # Filter to find the X_center and Y_center
    mask = np.zeros((M, N))
    mask[y1:y2, x1:x2] = 1
    filter = binary_FT * mask
    #imageShow(filter, 'FT Filter')


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
    imageShow(binary_phase, 'Binarized phase')


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
            #imageShow(phaseCompensate, 'phaseCompensate')

            minVal = np.amin(phaseCompensate)
            maxVal = np.amax(phaseCompensate)
            phase_sca = (phaseCompensate - minVal) / (maxVal - minVal)
            binary_phase = (phase_sca > 0.2)
            #imageShow(binary_phase, 'phaseCompensate')

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





'''
Auxiliary functions
'''

# Function to display an image
# Inputs:
# inp - The input complex field
# title - The title of the displayed image
def imageShow(inp, title):
    plt.imshow(inp, cmap='gray'), plt.title(title)  # image in gray scale
    plt.show()  # show image

    return

# Function to calculate the intensity representation of a given complex field
def intensity(inp, log):
    # Inputs:
    # inp - The input complex field
    # log - boolean variable to determine if a log representation is applied

    out = np.abs(inp)
    out = out * out

    if log == True:
        out = 20 * np.log(out)

    return out


# Spatial filtering for the FRS and ERS algorithms
def spatialFiltering(inp, M, N, upper):
    field_spec = np.fft.fftshift(inp)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)

    pow_espec = intensity(field_spec, True)

    # Defining width and hieght of point of maxima in the ROI determination
    # 10% is estimated to be the proper distance to separated an in-line hologram from an off-axis one in the Spectrum
    M_search = int(M - (M * 0.1))
    N_search = int(N - (N * 0.1))
    # print (M_search)

    # Location and size of the ROI for the +1 or -1 diffraction order info

    # Location of the center of the ROI in the DHM hologram spectrum
    if (upper == True):
        # Find the location of the maxima point in the upper part of the DHM hologram spectrum
        indI = np.unravel_index(np.argmax(pow_espec[1:(int(M / 2) - M_search), 1:N], axis=None),
                                pow_espec[1:(int(M / 2) - M_search), 1:N].shape)
        # print (indI)
    else:
        # Find the location of the maxima point in the lower part of the DHM hologram spectrum
        indI = np.unravel_index(np.argmax(pow_espec[1:M, 1:(int(N / 2) - N_search)], axis=None),
                                pow_espec[1:M, 1:(int(N / 2) - N_search)].shape)
        # print (indI)

    # Determination of the ROI size. To do this, we use the theoretical size of the +1 or -1 diffraction orders according...
    # ... to the distance of their centers to the DC coordinates (middle point in the DHM hologram spectrum).
    d = np.sqrt(np.power(indI[1] - M / 2, 2) + np.power(indI[0] - N / 2, 2))
    radius = d / 3

    # Mask building
    mask = np.zeros((M, N))
    for j in range(M):
        for i in range(N):
            if np.power(j - indI[1], 2) + np.power(i - indI[0], 2) < np.power(radius, 2):
                mask[i, j] = 1

    tmp = field_spec * mask
    # field_spec.fill(0)

    # pow_espec = intensity(tmp, True)
    # out = pow_espec

    # Coming back to spatial domain (retreiving filtered hologram)
    out = np.fft.ifftshift(tmp)
    out = np.fft.ifft2(out)
    out = np.fft.ifftshift(out)

    return out, indI, tmp


# Spatial filtering process for FCF implementation
def spatialFilteringCF(inp, M, N):
    field_spec = np.fft.fft2(inp)
    field_spec = np.fft.fftshift(field_spec)

    pow_espec = intensity(field_spec, True)

    # Finding the max peaks for +1 order in I or II quadrant
    height_half = round(int(M / 2))
    mask = np.zeros((M, N))
    mask[0:height_half - 1, 0:M] = 1
    field_spec_tem = pow_espec * mask
    maximum = np.amax(field_spec_tem)
    fy_max, fx_max = np.where(field_spec_tem == maximum)
    #print(fx_max, fy_max )

    # Determination of the ROI size. To do this, we use the theoretical size of the +1 or -1 diffraction orders according...
    # ... to the distance of their centers to the DC coordinates (middle point in the DHM hologram spectrum).
    d = np.sqrt(np.power(fx_max - M / 2, 2) + np.power(fy_max - N / 2, 2))
    radius = d / 3

    # Mask building
    mask = np.zeros((M, N), dtype=int)
    for p in range(0, N):
        for q in range(0, M):
            if np.sqrt((p - fy_max) ** 2 + (q - fx_max) ** 2) < radius:
                mask[p, q] = 1

    # Filtering the hologram
    tmp = field_spec * mask

    # Coming back to spatial domain (retrieving filtered hologram)
    out = np.fft.ifftshift(tmp)
    out = np.fft.ifft2(out)

    return out, fx_max, fy_max


# cost function for the CFS implementation
def costFunction(seeds, M, N, holo_filter, wavelength, dx, dy, X, Y, fx_0, fy_0, k):
    J = 0;
    theta_x = math.asin((fx_0 - seeds[0]) * wavelength / (M * dx))
    theta_y = math.asin((fy_0 - seeds[1]) * wavelength / (N * dy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    phase = compensation(holo_filter, ref_wave)

    phase = phase + math.pi
    sumIB = phase.sum()
    J = (M * N) - sumIB

    return J


# function for the compensation
def compensation(holo_filter, ref_wave):
    reconstruction = holo_filter * ref_wave
    phase = np.angle(reconstruction)

    return phase
    
