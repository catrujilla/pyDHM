# libraries
import numpy as np
import math

# phaseShifting package information
def info():
    print(
        "The available phase shifting algorithms are slightly off-axis synthetic reference (SOSR), "
        "blind shifting methods (BS2 and BS3) and the traditional full in-line 3-step, 4-step and 5-step algorithms")

    return


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
