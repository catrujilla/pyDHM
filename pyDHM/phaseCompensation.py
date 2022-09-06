"""
Title-->            Phase compensators package
Author-->           Carlos Trujillo, Ana Doblas and Raul Castaneda
Date-->             21/06/2021
Last modified-->    05/09/2022
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
import scipy
import sys
import cv2

def FRS(inp, upper, wavelength, dx, dy, s=2, step=10):
    '''
    # Function to compensate phase maps of off-axis DHM via the full ROI search algorithm.
    # Inputs:
    # inp - The input intensity (captured) hologram
    # upper - Boolean variable defining whether the ROI will be that on the upper part of the DHM hologram spectrum.
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    # s = 2 and steps = 10
    '''
    
    if step < s:
        print('Please, Enter a s value smaller than step')
        sys.exit()

    # determine if the hologram is off-axis
    orders, thresh = regime(inp)
    if orders < 3:
        print('FRS require an off-axis hologram')
        sys.exit()

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
    if step > s:
        print('Please, Enter a step value smaller than s')
        sys.exit()

    # determine if the hologram is off-axis
    orders, thresh = regime(inp)
    if orders < 3:
        print('ERS require an off-axis hologram')
        sys.exit()

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
    
    if scipy.__version__ == None:
        print('Please, install scipy.optimize library and import minimize function')
        sys.exit()

    # determine if the hologram is off-axis
    orders, thresh = regime(inp)
    if orders < 3:
        print('CFS require an off-axis hologram')
        sys.exit()

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
    step = 1
    res = minimize(costFunction, seeds, args=(M, N, holo_filter, wavelength, dx, dy, X, Y, fx_0, fy_0, k),
                   method='TNC', bounds=((seeds[0] - step, seeds[0] + step), (seeds[1] - step, seeds[1] + step)), tol=1e-3)
    print("Minimization process finished.")
    fx_max = res.x[0]
    fy_max = res.x[1]
    print('x: ',fx_max)
    print('y: ',fy_max)
    
    # Best phase compensation
    print("Phase compensation started....")
    theta_x = math.asin((fx_0 - fx_max) * wavelength / (M * dx))
    theta_y = math.asin((fy_0 - fy_max) * wavelength / (N * dy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    comp_phase = holo_filter * ref_wave

    print("Phase compensation finished.")

    return comp_phase


def CNT(inp, wavelength, dx, dy, x1=None, x2=None, y1=None, y2=None, spatialFilter=None):
    '''
    # Function to compensate phase maps of image plane off-axis DHM, operating in non-telecentric regimen
    # Inputs:
    # inp - The input intensity (captured) hologram
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    # x1 - Coordinate x1 for the rectangle (upper left corner)
    # y1 - Coordinate y1 for the rectangle (upper left corner)
    # x2 - Coordinate x2 for the rectangle (lower right corner)
    # y2 - Coordinate y2 for the rectangle (lower right corner)
    # spatialFilter - The approach to compute the spatial filter, two options available      and sfmr
    '''

    wavelength = wavelength * 0.000001
    dx = dx * 0.000001
    dy = dy * 0.000001

    # Retrieving the input shape
    inp = np.array(inp)
    M, N = inp.shape
    k = (2 * math.pi) / wavelength

    # Creating a mesh-grid to operate in world coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY

    # The spatial filtering process is executed
    print("Spatial filtering process started.....")
    if x1 is None and x2 is None and y1 is None and y2 is None:
        if spatialFilter == 'sfmr':
            Xcenter, Ycenter, holo_filter, ROI_array = spatialFilterinCNT(inp, M, N)
        else:
            print("Please, indicate as option for the spatialFilter: 'sfmr' ")
            sys.exit()
    else:
        if spatialFilter == 'sfr':
            Xcenter, Ycenter, holo_filter, ROI_array = spatialFilterinCNT_II(inp, M, N, x1, y1, x2, y2)
        else:
            print("Please, indicate as option for the spatialFilter: 'sfr' or introduce the rectangle coordinates")
            sys.exit()
    print("Spatial filtering process finished.")

    # Fourier transform to the hologram filtered
    ft_holo = FT(holo_filter)
    FT_display = intensity(ft_holo, False)
    #imageShow(FT_display, 'FT Filtered')

    # reference wave for the first compensation (global linear compensation)
    ThetaXM = math.asin((N / 2 - Xcenter) * wavelength / (M * dx))
    ThetaYM = math.asin((M / 2 - Ycenter) * wavelength / (N * dy))
    reference = np.exp(1j * k * (math.sin(ThetaXM) * X * dx + math.sin(ThetaYM) * Y * dy))

    # First compensation
    comp_phase = holo_filter * reference
    phase_c = phase(comp_phase)

    # show the first compensation
    minVal = np.amin(phase_c)
    maxVal = np.amax(phase_c)
    phase_normalized = (phase_c - minVal) / (maxVal - minVal)
    binary_phase = (phase_normalized > 0.2)
    imageShow(binary_phase, 'Binarized phase')


    # creating the new reference wave to eliminate the circular phase factors
    m = abs(ROI_array[2] - ROI_array[0])
    n = abs(ROI_array[3] - ROI_array[1])
    Cx = np.power((M * dx), 2)/(wavelength * m)
    Cy = np.power((N * dy), 2)/(wavelength * n)
    cur = (Cx + Cy)/2

    print("Carefully determine the center of the circular phase factor in the Binarized Image...")
    p = input("Enter the pixel position X_cent of the center of circular phase map on x axis ")
    q = input("Enter the pixel position Y_cent of the center of circular phase map on y axis ")
    f = ((M/2) - int(p))/2
    g = ((N/2) - int(q))/2
    print("Phase compensation started....")

    cont = 0
    sum_max = 0
    s = 100
    step = 50
    perc = 40/100

    arrayCurvature = np.arange(cur - (cur*perc), cur + (cur*perc), perc/6)
    arrayXcenter = np.arange(f - s, f + s, step)
    arrayYcenter = np.arange(g - s, g + s, step)
    for curTemp in arrayCurvature:
        for fTemp in arrayXcenter:
            for gTemp in arrayYcenter:
                cont = cont + 1
                phi_spherical = (np.power(X - fTemp, 2) * np.power(dx, 2) / curTemp) + (
                np.power(Y - gTemp, 2) * np.power(dy, 2) / curTemp)
                phi_spherical = math.pi * phi_spherical / wavelength
                phi_spherical = np.exp(-1j * phi_spherical)

                phaseCompensate = comp_phase * phi_spherical
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
                    f_out = fTemp
                    g_out = gTemp
                    cur_out = curTemp
                    sum_max = sum

    #print("after first search ", f_out, g_out, cur_out)

    cont = 0
    sum_max = 0
    s = 10
    step = 2
    perc = 0.1
    arrayXcenter = np.arange(f_out - s, f_out + s, step)
    arrayYcenter = np.arange(g_out - s, g_out + s, step)
    arrayCurvature = np.arange(cur_out - (cur_out*perc), cur_out + (cur_out*perc), 0.01)
    #arrayCurvature = np.arange(1.003, 1.03, 0.01)

    for curTemp in arrayCurvature:
        for fTemp in arrayXcenter:
            for gTemp in arrayYcenter:
                #print(curTemp)
                cont = cont + 1
                phi_spherical = (np.power(X - fTemp, 2) * np.power(dx, 2) / curTemp) + (
                    np.power(Y - gTemp, 2) * np.power(dy, 2) / curTemp)
                phi_spherical = math.pi * phi_spherical / wavelength
                phi_spherical = np.exp(-1j * phi_spherical)

                phaseCompensate = comp_phase * phi_spherical
                phaseCompensate = np.angle(phaseCompensate)
                #imageShow(phaseCompensate, 'phaseCompensate')

                minVal = np.amin(phaseCompensate)
                maxVal = np.amax(phaseCompensate)
                phase_sca = (phaseCompensate - minVal) / (maxVal - minVal)
                binary_phase = (phase_sca > 0.2)
                #imageShow(binary_phase, 'phaseCompensate')

                # Applying the summation and thresholding metric
                sum = np.sum(np.sum(binary_phase))
                #print(sum, curTemp)
                if (sum > sum_max):
                    f_out = fTemp
                    g_out = gTemp
                    cur_out = curTemp
                    sum_max = sum

    phi_spherical = (np.power(X - f_out, 2) * np.power(dx, 2) / cur_out) + (
            np.power(Y - g_out, 2) * np.power(dy, 2) / cur_out)
    phi_spherical = math.pi * phi_spherical / wavelength
    phi_spherical = np.exp(-1j * phi_spherical)
    phaseCompensate = comp_phase * phi_spherical

    
    #print("after fine compensation", f_out, g_out, cur_out)

    print("Phase compensation finished.")

    return phaseCompensate

'''
Auxiliary functions
'''

# Function to determine if the holograms is off-axis or not
def regime(inp):
    holoFT = np.float32(inp)
    fft_holo = cv2.dft(holoFT, flags=cv2.DFT_COMPLEX_OUTPUT)
    fft_holo = np.fft.fftshift(fft_holo)
    fft_holo_image = 20 * np.log(cv2.magnitude(fft_holo[:, :, 0], fft_holo[:, :, 1]))
    minVal = np.amin(np.abs(fft_holo_image))
    maxVal = np.amax(np.abs(fft_holo_image))
    fft_holo_image = cv2.convertScaleAbs(fft_holo_image, alpha=255.0 / (maxVal - minVal),
                                         beta=-minVal * 255.0 / (maxVal - minVal))

    # cv2.imshow('Binary image_resize', fft_holo_image)
    # cv2.waitKey(0)

    # apply binary thresholding
    ret, thresh = cv2.threshold(fft_holo_image, 200, 255, cv2.THRESH_BINARY)
    #cv2.imshow('Binary image', thresh)
    thresh_rize = cv2.resize(thresh, (1024, 1024))
    #cv2.imshow('Binary image_resize', thresh_rize)
    #cv2.waitKey(0)


    contours, hierarchy = cv2.findContours(image=thresh_rize, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_SIMPLE)
    # draw contours on the original image
    fft_holo_image = cv2.resize(thresh, (1024, 1024))
    image_copy = fft_holo_image.copy()
    cv2.drawContours(image=image_copy, contours=contours, contourIdx=-1, color=(0, 255, 0), thickness=2,
                     lineType=cv2.LINE_AA)
    # cv2.imshow('None approximation', image_copy)
    # cv2.waitKey(0)
    orders = len(contours)
    # print(orders)
    return orders, thresh

# Spatial filtering process for rectangular selection for CNT
def spatialFilterinCNT_II(inp, M, N, x1, y1, x2, y2):
    ROI_array = np.zeros(4)

    # Fourier transform of the hologram
    FT = np.fft.fft2(inp)
    FT = np.fft.fftshift(FT)

    # Filter to find the X_center and Y_center
    mask = np.zeros((M, N))
    mask[y1:y2, x1:x2] = 1

    # computing the center of the rectangle mask (X_center and Y_center)
    Xcenter = x1 + (x2 - x1)/2
    Ycenter = y1 + (y2 - y1)/2

    # spatial filter
    filter = FT * mask
    holo_filter = np.fft.ifftshift(filter)
    holo_filter = np.fft.ifft2(holo_filter)

    # creating array for coordinates x1, y1, x2, y2
    ROI_array[0] = x1
    ROI_array[1] = y1
    ROI_array[2] = x1 + x2
    ROI_array[3] = y1 + y2

    return Xcenter, Ycenter, holo_filter, ROI_array

# Spatial filtering process - manual selection for CNT
def spatialFilterinCNT(inp, M, N):
    ROI_array = np.zeros(4)
    holoFT = np.float32(inp)  # convertion of data to float
    fft_holo = cv2.dft(holoFT, flags=cv2.DFT_COMPLEX_OUTPUT)  # FFT of hologram
    fft_holo = np.fft.fftshift(fft_holo)
    fft_holo_image = 20 * np.log(cv2.magnitude(fft_holo[:, :, 0], fft_holo[:, :, 1]))  # logaritm scale FFT
    minVal = np.amin(np.abs(fft_holo_image))
    maxVal = np.amax(np.abs(fft_holo_image))
    fft_holo_image = cv2.convertScaleAbs(fft_holo_image, alpha=255.0 / (maxVal - minVal),
                                         beta=-minVal * 255.0 / (maxVal - minVal))

    ROI = cv2.selectROI(fft_holo_image, fromCenter=True)  # module to  ROI
    # imCrop = fft_holo_image[int(ROI[1]):int(ROI[1] + ROI[3]), int(ROI[0]):int(ROI[0] + ROI[2])]
    x1_ROI = int(ROI[1])
    y1_ROI = int(ROI[0])
    x2_ROI = int(ROI[1] + ROI[3])
    y2_ROI = int(ROI[0] + ROI[2])
    ROI_array[0] = x1_ROI
    ROI_array[1] = y1_ROI
    ROI_array[2] = x2_ROI
    ROI_array[3] = y2_ROI

    # computing the center of the rectangle mask
    Ycenter = x1_ROI + (x2_ROI - x1_ROI)/2
    Xcenter = y1_ROI + (y2_ROI - y1_ROI)/2

    holo_filter = np.zeros((M, N, 2))
    holo_filter[x1_ROI:x2_ROI, y1_ROI: y2_ROI] = 1
    holo_filter = holo_filter * fft_holo
    holo_filter = np.fft.ifftshift(holo_filter)
    holo_filter = cv2.idft(holo_filter, flags=cv2.DFT_INVERSE)

    holo_filter_real = holo_filter[:, :, 0]
    holo_filter_imag = holo_filter[:, :, 1]
    holo_filter = np.zeros((M, N), complex)
    for p in range(M):
        for q in range(N):
            holo_filter[p, q] = complex(holo_filter_real[p, q], holo_filter_imag[p, q])

    return Xcenter, Ycenter, holo_filter, ROI_array

# Function to display an image
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

# Function to calcule the amplitude representation of a given complex field
def amplitude(inp, log):
    out = np.abs(inp)

    if log == True:
        out = 20 * np.log(out)

    return out

# Function to calcule the phase representation of a given complex field using the
def phase(inp):
    out = np.angle(inp)

    return out
    
# Function to calcule the Fourier transform of a given field using the
def FT(inp):
    FT = np.fft.fft2(inp)
    FT = np.fft.fftshift(FT)

    return FT

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
    
