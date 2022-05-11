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
import phaseShifting._init_ as phs
from scipy.optimize import minimize


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
    pow_espec = phs.intensity(field_spec, True)

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

    print("Synthetic building of the reference wave started....")

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
    print("Synthetic reference wave built.")

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

    print("Phase compensation finished.")

    return comp_phase


def BPS2(Inp0, Inp1, wavelength, dx, dy):
    '''
    # Function to recover the phase information of a sample from two slightly off-axis acquisitions.
    # Inputs:
    # InpX - The input intensities
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dx, dy- Pixel dimensions of the camera sensor used for recording the holograms
    '''

    # retrieving the input shape
    inp = np.array(Inp0)
    M, N = inp.shape

    # FT hologram
    field_spec = np.fft.fft2(inp)
    field_spec = np.fft.fftshift(field_spec)
    pow_espec = phs.intensity(field_spec, True)

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
    res = minimize(phs.costFunction2, test_theta, args=(Inp0, Inp1, fx_max_d1, fy_max_d1, fx_max_d2,
                                                        fy_max_d2), method='Cobyla', tol=1e-6)
    print("Minimization process finished.")

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
    out, fx_max, fy_max = phs.spatialFiltering(FTd2, M, N)

    # compensation process
    print("Phase compensation started....")
    # computing the digital reference wave
    ref_wave = phs.referenceWave(M, N, wavelength, dx, dy, fx_max, fy_max)
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

    # Retrieving the input shape
    inp = np.array(Inp0)
    M, N = inp.shape

    # FT hologram
    field_spec = np.fft.fft2(inp)
    field_spec = np.fft.fftshift(field_spec)
    pow_espec = phs.intensity(field_spec, True)

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
    res = minimize(phs.costFunction3, test_theta, args=(Inp0, Inp1, Inp2, fx_max_d1, fy_max_d1, fx_max_d2,
                                                    fy_max_d2), method='Cobyla', tol=1e-6)
    print("Minimization process finished.")

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
    out, fx_max, fy_max = phs.spatialFiltering(FTd, M, N)

    # compensation process
    print("Phase compensation started....")
    # computing the digital reference wave
    ref_wave = phs.referenceWave(M, N, wavelength, dx, dy, fx_max, fy_max)
    comp_phase = out * ref_wave
    print("Phase compensation finished.")

    return comp_phase


def PS5(Inp0, Inp1, Inp2, Inp3, Inp4):
    '''
    # Function to recover the phase information of a sample from five DHM in-axis acquisitions holograms
    # Inputs:
    # inpX - The input intensities (captured) pi/2 phase-shifted holograms
    '''

    # Retrieving the input shape
    inp0 = np.array(Inp0)
    inp1 = np.array(Inp1)
    inp2 = np.array(Inp2)
    inp3 = np.array(Inp3)
    inp4 = np.array(Inp4)

    print("Phase compensation started....")
    # computing the compensation
    comp_phase = np.arctan((2*(inp3-inp1))/(2*inp2-inp0-inp4))
    print("Phase compensation finished.")

    return comp_phase


def PS4(Inp0, Inp1, Inp2, Inp3):
    '''
    # Function to recover the phase information of a sample from four DHM in-axis acquisitions holograms
    # Inputs:
    # inpX - The input intensities (captured) pi/2 phase-shifted holograms
    '''

    # Retrieving the input shape
    inp0 = np.array(Inp0)
    inp1 = np.array(Inp1)
    inp2 = np.array(Inp2)
    inp3 = np.array(Inp3)


    # compensation process
    print("Phase compensation started....")
    # computing the compensation
    comp_phase = np.arctan((inp3-inp1)/(inp2-inp0))
    print("Phase compensation finished.")

    return comp_phase


def PS3(Inp0, Inp1, Inp2):
    '''
    # Function to recover the phase information of a sample from three DHM in-axis acquisitions holograms
    # Inputs:
    # inpX - The input intensities (captured) pi/3 phase-shifted holograms
    '''

    # Retrieving the input shape
    inp0 = np.array(Inp0)
    inp1 = np.array(Inp1)
    inp2 = np.array(Inp2)


    print("Phase compensation started....")

    # computing the compensation
    comp_phase = np.arctan((np.sqrt(3)*(inp2-inp0))/((inp2+inp0)-(2*inp1)))
    print("Phase compensation finished.")

    return comp_phase