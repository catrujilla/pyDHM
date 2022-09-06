"""
Title-->            Phase phaseShifting package
Author-->           Ana Doblas, Carlos Trujillo, Raul Castaneda,
Date-->             21/02/2021
Last modified-->    11/04/2022
Groups-->           University of Memphis -  Optical Imaging Research laboratory (OIRL)
                    EAFIT University - Applied Optics Group
Abstract -->        Script for implementing six phase shifting algorithms implemented in pyDHM:
                    SOSR, BPS2, BPS3, PS5, PS4 and PS3.
Links-->            More information can be found on the website
                    https://github.com/catrujilla/pyDHM
"""

# libraries
import numpy as np
import math
from math import pi
from scipy.optimize import minimize
import scipy
import sys
import cv2

def SOSR(inp0, inp1, inp2, inp3, upper, wavelength, dx, dy, s=1, steps=4):
    '''
    # Function to recover the phase information of a sample from 4 DHM slightly off-axis acquisitions (quadrature Phase Shifting method).
    # The idea was first proposed by Ferraro, we upgraded the method to accurately calculate the synthetic reference wave
    # to obtain a fully compensated phase maps.
    # Inputs:
    # inpX - The input intensities (captured) pi/2 phase-shifted holograms
    # upper - Boolean variable defining whether the brightest point will be that on the upper part of the DHM hologram spectrum.
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the holograms
    '''

    if steps < s:
        print('Please, Enter a s value smaller than step')
        sys.exit()

    # Retrieving the input shape

    inp0 = inp0 - np.average(inp0)
    inp1 = inp1 - np.average(inp1)
    inp2 = inp2 - np.average(inp2)
    inp3 = inp3 - np.average(inp3)
    M, N = inp0.shape

    # Creating a mesh_grid to operate in world-coordinates
    x = np.arange(0, N, 1)  # array x
    y = np.arange(0, M, 1)  # array y
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')  # meshgrid XY

    # Phase Shifting process

    # Location of the brightest point in the first phase-shifted DHM hologram spectrum
    # First, the need to Fourier transform and biuld the power spectrum
    field_spec = np.fft.fftshift(inp0)
    field_spec = np.fft.fft2(field_spec)
    field_spec = np.fft.fftshift(field_spec)
    pow_espec = intensity(field_spec, True)

    # Defining width and height of point of maxima in the ROI determination
    # 10% is estimated to be the proper distance to separated an in-line hologram from an off-axis one in the Spectrum
    M_search = int(M - (M * 0.1))
    N_search = int(N - (N * 0.1))
    # print (M_search)

    # Determining the point of maxima in the power spectrum
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

    # DC coordinates (centers of the DHM hologram Spectrums)
    fx_0, fy_0 = N / 2, M / 2

    # Coordinates of maxima points in the DHM hologram spectrums
    fx_1, fy_1 = indI[1] + 1, indI[0] + 1  # +1 is the correction of the offset added by numpy?

    # print (fx_0, fy_0, fx_1, fy_1)

    sum_max = 0
    arrayX = np.linspace(fx_1 - s, fx_1 + s, steps)
    arrayY = np.linspace(fy_1 - s, fy_1 + s, steps)

    k = (2 * pi) / wavelength

    print("Synthetic building of the reference wavefront started....")

    for fx_temp in arrayX:
        for fy_temp in arrayY:
            # print (fx_temp, fy_temp)
            theta_x = math.asin((fx_0 - fx_temp) * wavelength / (N * dx))
            theta_y = math.asin((fy_0 - fy_temp) * wavelength / (M * dy))

            R0 = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
            R1 = 1j * np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
            R2 = (-1) * np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
            R3 = (-1j) * np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))

            # Four-Quarter Phase Shifting (Ferraro's proposal)
            Isynth = (R0 * inp0) + (R1 * inp1) + (R2 * inp2) + (R3 * inp3)

            # module 2-pi phase retrieval
            phase = np.angle(Isynth)

            # Thresholding process
            minVal = np.amin(phase)
            maxVal = np.amax(phase)
            phase_sca = (phase - minVal) / (maxVal - minVal)
            binary_phase = (phase_sca < 0.1)

            # plt.imshow(phase, cmap='gray'), plt.title("phase")  # image in gray scale
            # plt.show()  # show image

            # Applying the summation and thresholding metric
            sum = np.sum(np.sum(binary_phase))
            if sum > sum_max:
                x_max_out = fx_temp;
                y_max_out = fy_temp;
                sum_max = sum;
            # print (sum)

    # print(x_max_out, y_max_out)
    print("Synthetic reference wavefront built.")

    # Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * dx))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * dy))
    R0 = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    R1 = 1j * np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    R2 = (-1) * np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    R3 = (-1j) * np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))

    # Four-Quarter Phase Shifting (Ferraro's proposal)
    Isynth = (R0 * inp0) + (R1 * inp1) + (R2 * inp2) + (R3 * inp3)

    comp_phase = Isynth

    print("Phase-shifting reconstruction finished.")

    return comp_phase


def BPS2(Inp0, Inp1, wavelength, dx, dy):
    '''
    # Function to recover the phase information of a sample from two slightly off-axis acquisitions.
    # Inputs:
    # InpX - The input intensities
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy- Pixel dimensions of the camera sensor used for recording the holograms
    '''

    if scipy.__version__ == None:
        print('Please, install scipy library and import minimize function')
        sys.exit()

    # retrieving the input shape
    inp = np.array(Inp0)
    M, N = inp.shape

    # FT hologram
    field_spec = np.fft.fft2(inp)
    field_spec = np.fft.fftshift(field_spec)
    pow_espec = intensity(field_spec, True)

    # finding the max peaks for +1 order in I or II quadrant
    height_half = round(int(M / 2))
    mask = np.zeros((M, N))
    mask[0:height_half - 1, 0:N] = 1
    field_spec_tem = pow_espec * mask
    maximum = np.amax(field_spec_tem)
    fy_max_d1, fx_max_d1 = np.where(field_spec_tem == maximum)

    # finding the max peaks for -1 order in III or IV quadrant
    height_half = round(int(M / 2))
    mask = np.zeros((M, N))
    mask[height_half + 1:M, 0:N] = 1
    field_spec_tem = pow_espec * mask
    maximum = np.amax(field_spec_tem)
    fy_max_d2, fx_max_d2 = np.where(field_spec_tem == maximum)

    # creating seeds for the minimization
    test_theta = np.random.randint(0, 360, (1, 2))
    test_theta = test_theta * math.pi / 180
    test_theta = [test_theta[0, 0], test_theta[0, 1]]

    # minimization
    print("Minimization process started.....")
    res = minimize(costFunction2, test_theta, args=(Inp0, Inp1, fx_max_d1, fy_max_d1, fx_max_d2,
                                                        fy_max_d2), method='Cobyla', tol=1e-6)
    print("Minimization process finished.")

    print("Phase compensation started....")
    # computing d1 and d2
    x = res.x
    theta1 = x[0]
    theta2 = x[1]
    ones = np.ones((2, 2))
    Matrix = [[1, np.exp(1j*theta1)],
         [1, np.exp(1j*theta2)]]
    Matrix = Matrix * ones
    Minv = np.linalg.inv(Matrix)
    d2 = Minv[1, 0]*Inp0 + Minv[1, 1]*Inp1

    # FFT d2
    FTd2 = np.fft.fft2(d2)
    FTd2 = np.fft.fftshift(FTd2)

    # spatial filtering
    out, fx_max, fy_max = spatialFiltering(FTd2, M, N)

    # compensation process
    # computing the digital reference wave
    ref_wave = referenceWave(M, N, wavelength, dx, dy, fx_max, fy_max)
    comp_phase = out * ref_wave
    # introducir metodo de compensación
    print("Phase compensation finished.")

    return comp_phase


def BPS3(Inp0, Inp1, Inp2, wavelength, dx, dy):
    '''
    # Function to recover the phase information of a sample from three slightly off-axis acquisitions.
    # Inputs:
    # InpX - The input intensities (captured)
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy - Pixel dimensions of the camera sensor used for recording the hologram
    '''
    
    if scipy.__version__ == None:
        print('Please, install scipy library and import minimize function')
        sys.exit()


    # Retrieving the input shape
    inp = np.array(Inp0)
    M, N = inp.shape

    # FT hologram
    field_spec = np.fft.fft2(inp)
    field_spec = np.fft.fftshift(field_spec)
    pow_espec = intensity(field_spec, True)

    # finding the max peaks for +1 order in I or II quadrant
    height_half = round(int(M / 2))
    mask = np.zeros((M, N))
    mask[0:height_half - 1, 0:N] = 1
    field_spec_tem = pow_espec * mask
    maximum = np.amax(field_spec_tem)
    fy_max_d1, fx_max_d1 = np.where(field_spec_tem == maximum)

    # finding the max peaks for -1 order in III or IV quadrant
    height_half = round(int(M / 2))
    mask = np.zeros((M, N))
    mask[height_half + 1:M, 0:N] = 1
    field_spec_tem = pow_espec * mask
    maximum = np.amax(field_spec_tem)
    fy_max_d2, fx_max_d2 = np.where(field_spec_tem == maximum)

    # creating seeds for the minimization
    test_theta = np.random.randint(0, 360, (1, 3))
    test_theta = test_theta * math.pi / 180
    test_theta = [test_theta[0, 0], test_theta[0, 1], test_theta[0, 2]]

    # minimization
    print("Minimization process started.....")
    res = minimize(costFunction3, test_theta, args=(Inp0, Inp1, Inp2, fx_max_d1, fy_max_d1, fx_max_d2,
                                                    fy_max_d2), method='Cobyla', tol=1e-6)
    print("Minimization process finished.")

    print("Phase compensation started....")
    # computing d1, d2 and d3
    x = res.x
    theta1 = x[0]
    theta2 = x[1]
    theta3 = x[2]
    ones = np.ones((3, 3))
    Matrix = [[1, np.exp(1j * theta1), np.exp(-1j * theta1)],
         [1, np.exp(1j * theta2), np.exp(-1j * theta2)],
         [1, np.exp(1j * theta3), np.exp(-1j * theta3)]]
    Matrix = Matrix * ones
    Minv = np.linalg.inv(Matrix)
    d = Minv[1, 0] * Inp0 + Minv[1, 1] * Inp1 + Minv[1, 2] * Inp2


    # FFT d
    FTd = np.fft.fft2(d)
    FTd = np.fft.fftshift(FTd)

    # spatial filtering
    out, fx_max, fy_max = spatialFiltering(FTd, M, N)

    # compensation process
    # computing the digital reference wave
    ref_wave = referenceWave(M, N, wavelength, dx, dy, fx_max, fy_max)
    comp_phase = out * ref_wave
    print("Phase compensation finished.")

    return comp_phase


def PS5(Inp0, Inp1, Inp2, Inp3, Inp4):
    '''
    # Function to recover the phase information of a sample from five DHM in-axis acquisitions holograms
    # Inputs:
    # inpX - The input intensities (captured) pi/2 phase-shifted holograms
    '''
    
    # determine if the hologram is on-axis
    ret, thresh = regime(Inp0)
    if ret != 1:
        print('PS5 require on-axis holograms')
        sys.exit()

    # Retrieving the input shape
    inp0 = np.array(Inp0)
    inp1 = np.array(Inp1)
    inp2 = np.array(Inp2)
    inp3 = np.array(Inp3)
    inp4 = np.array(Inp4)

    print("Phase-shifing reconstruction started....")
    # computing the compensation
    comp_phase = np.arctan((2*(inp3-inp1))/(2*inp2-inp0-inp4))
    print("Phase-shifing reconstruction finished.")

    return comp_phase


def PS4(Inp0, Inp1, Inp2, Inp3):
    '''
    # Function to recover the phase information of a sample from four DHM in-axis acquisitions holograms
    # Inputs:
    # inpX - The input intensities (captured) pi/2 phase-shifted holograms
    '''
    
    # determine if the hologram is on-axis
    ret, thresh = regime(Inp0)
    if ret != 1:
        print('PS4 require on-axis holograms')
        sys.exit()


    # Retrieving the input shape
    inp0 = np.array(Inp0)
    inp1 = np.array(Inp1)
    inp2 = np.array(Inp2)
    inp3 = np.array(Inp3)


    # compensation process
    print("Phase-shifing reconstruction started....")
    # computing the compensation
    comp_phase = np.arctan((inp3-inp1)/(inp2-inp0))
    print("Phase-shifing reconstruction finished.")

    return comp_phase


def PS3(Inp0, Inp1, Inp2):
    '''
    # Function to recover the phase information of a sample from three DHM in-axis acquisitions holograms
    # Inputs:
    # inpX - The input intensities (captured) pi/3 phase-shifted holograms
    '''

    # determine if the hologram is on-axis
    ret, thresh = regime(Inp0)
    if ret != 1:
        print('PS3 require on-axis holograms')
        sys.exit()

    # Retrieving the input shape
    inp0 = np.array(Inp0)
    inp1 = np.array(Inp1)
    inp2 = np.array(Inp2)


    print("Phase-shifing reconstruction started....")

    # computing the compensation
    comp_phase = np.arctan((np.sqrt(3)*(inp2-inp0))/((inp2+inp0)-(2*inp1)))
    print("Phase-shifing reconstruction finished.")

    return comp_phase
    
    
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

# Function to calculate the intensity representation of a given complex field
def intensity(inp, log):
    # Function to calcule the intensity representation of a given complex field
    # Inputs:
    # inp - The input complex field
    # log - boolean variable to determine if a log representation is applied
    out = np.abs(inp)
    out = out * out

    if log == True:
        out = 20 * np.log(out)

    return out


# cost function for the BPS2 implementation
def costFunction2(test_theta, hologram1, hologram2, fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2):
    theta1 = test_theta[0]
    theta2 = test_theta[1]
    ones = np.ones((2, 2))
    Matrix = [[1, np.exp(1j * theta1)],
         [1, np.exp(1j * theta2)]]
    Matrix = Matrix * ones
    Minv = np.linalg.inv(Matrix)
    d1 = Minv[0, 0] * hologram1 + Minv[0, 1] * hologram2

    # FFT d1
    FTd1 = np.fft.fft2(d1)
    FTd1 = np.fft.fftshift(FTd1)

    # calculate J
    J = 1 - (np.abs(FTd1[fx_max_d1, fy_max_d1]) - np.abs(FTd1[fx_max_d2, fy_max_d2])) / \
        (np.abs(FTd1[fx_max_d2, fy_max_d2]) + np.abs(FTd1[fx_max_d1, fy_max_d1]))
    J = J[0]
    return J


# cost function for the BPS3 implementation
def costFunction3(test_theta, hologram1, hologram2, hologram3, fx_max_d1, fy_max_d1, fx_max_d2, fy_max_d2):
    theta1 = test_theta[0]
    theta2 = test_theta[1]
    theta3 = test_theta[2]
    ones = np.ones((3, 3))
    Matrix = [[1, np.exp(1j*theta1), np.exp(-1j*theta1)],
         [1, np.exp(1j*theta2), np.exp(-1j*theta2)],
         [1, np.exp(1j*theta3), np.exp(-1j*theta3)]]
    Matrix = Matrix * ones
    Minv = np.linalg.inv(Matrix)
    d3 = Minv[2, 0]*hologram1 + Minv[2, 1]*hologram2 + Minv[2, 2]*hologram3

    # FFT d3
    FTd3 = np.fft.fft2(d3)
    FTd3 = np.fft.fftshift(FTd3)

    # calculate J
    J = 1 + (1 - (np.abs(FTd3[fx_max_d2, fy_max_d2]) - np.abs(FTd3[fx_max_d1, fy_max_d1])) /
             (np.abs(FTd3[fx_max_d2, fy_max_d2]) + np.abs(FTd3[fx_max_d1, fy_max_d1])))
    J = J[0]
    return J


# Spatial filtering process for blind phase-shifting methods
def spatialFiltering(inp, M, N):
    # Find the max peaks I or II quadrant
    height_half = round(int(M / 2))
    mask = np.zeros((M, N))
    mask[0:height_half - 1, 0:N] = 1
    field_spec_tem = inp * mask
    maximum = np.amax(field_spec_tem)
    fy_max, fx_max = np.where(field_spec_tem == maximum)

    # Determination of the ROI size. To do this, we use the theoretical size of the +1 or -1 diffraction orders according...
    # ... to the distance of their centers to the DC coordinates (middle point in the DHM hologram spectrum).
    # d = np.sqrt(np.power(fx_max - M / 2, 2) + np.power(fy_max - N / 2, 2))
    # radius = d / 3
    radius = 100

    # Mask building
    mask = np.zeros((M, N), dtype=int)
    for p in range(0, N):
        for q in range(0, M):
            if np.sqrt((p - fy_max) ** 2 + (q - fx_max) ** 2) < radius:
                mask[p, q] = 1

    # Filtering the hologram
    tmp = inp * mask

    # Coming back to spatial domain (retrieving filtered hologram)
    out = np.fft.ifftshift(tmp)
    out = np.fft.ifft2(out)

    return out, fx_max, fy_max


# function to create the reference wave
def referenceWave(M, N, wavelength, dx, dy, fx_max, fy_max):
    x = np.arange(0, M, 1)  # array x
    y = np.arange(0, N, 1)  # array y
    X, Y = np.meshgrid(x - (M / 2), y - (N / 2), indexing='xy')  # meshgrid XY
    fx_0 = M / 2
    fy_0 = N / 2
    k = (2 * math.pi) / wavelength
    theta_x = math.asin((fx_0 - fx_max) * wavelength / (M * dx))
    theta_y = math.asin((fy_0 - fy_max) * wavelength / (N * dy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dx) + (math.sin(theta_y) * Y * dy)))
    return ref_wave
