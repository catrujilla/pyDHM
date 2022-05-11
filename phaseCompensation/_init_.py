# libraries
import numpy as np
import math


# phaseCompensators package information
def info():
    print(
        "The available phase compensation algorithms are Full ROI Search (FRS), Efficient ROI Search (ERS) "
        "the Cost Function Search (CFS) and the compensation no-tele (CNT)")

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
