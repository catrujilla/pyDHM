

import numpy as np
import math
from scipy.optimize import minimize


def cut_hologram(hologram):
    height, width = hologram.shape
    if width > height:
        cut = round(int(width - height) / 2)
        hologram = hologram[0:height, cut:(width - cut)]
    elif width < height:
        cut = round(int(height - width) / 2)
        hologram = hologram[cut:(height - cut), 0:width]
    else:
        hologram = hologram
    height, width = hologram.shape
    return hologram, height, width


def spatial_filter(hologram, height, width, diameter):

    # FFT of the hologram
    fft = np.fft.fft2(hologram)
    fft = np.fft.fftshift(fft)
    fft_image = 20 * np.log(np.abs(fft))

    # Find the max peaks I or II quadrant
    height_half = round(int(height / 2))
    maskII = np.zeros((height, width))
    maskII[0:height_half - 1, 0:width] = 1
    fft_images_II = fft_image * maskII
    maximum = np.amax(fft_images_II)
    fy_max, fx_max = np.where(fft_images_II == maximum)

    # creating mask circular or square
    mask = np.zeros((height, width), dtype=int)
    for p in range(0, height):
        for q in range(0, width):
            if np.sqrt((p - fy_max) ** 2 + (q - fx_max) ** 2) < diameter:
                mask[p, q] = 1
    '''
    limit = 50
    x1 = round(int(fx_max - limit))
    x2 = round(int(fx_max + limit))
    y1 = round(int(fy_max - limit))
    y2 = round(int(fy_max + limit))
    mask[x1:x2, y1:y2] = 1
    '''

    # Filter hologram
    hologram_filter = fft * mask
    hologram_filter = np.fft.ifftshift(hologram_filter)
    hologram_filter = np.fft.ifft2(hologram_filter)

    return fx_max, fy_max, hologram_filter


def cost_function(peaks, height, width, hologram_filter, wave_length, deltaX, deltaY):
    J = 0;
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY
    fx_0 = height / 2
    fy_0 = width / 2
    k = (2 * math.pi) / wave_length
    theta_x = math.asin((fx_0 - peaks[0]) * wave_length / (height * deltaX))
    theta_y = math.asin((fy_0 - peaks[1]) * wave_length / (width * deltaY))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    phase = numerical_reconstructions(hologram_filter, ref_wave)
    '''
    phase = phase + math.pi
    minVal = np.amin(np.abs(phase))
    maxVal = np.amax(np.abs(phase))
    phase_image = cv2.convertScaleAbs(phase, alpha=255.0 / (maxVal - minVal),
                                         beta=-minVal * 255.0 / (maxVal - minVal))
    ret, ib = cv2.threshold(phase_image, 10, 1, cv2.THRESH_BINARY)
    # image_show(ib, 'Phase binarized')
    ibinary = np.ones((height, width), dtype=int)
    ibinary = ib * ibinary
    sumIB = ibinary.sum()
    J = (height * width) - sumIB
    '''
    phase = phase + math.pi
    sumIB = phase.sum()
    J = (height * width)*math.pi - sumIB
    #print('fx, fy: ', peaks[0], peaks[1])
    return J


def reference_wave(height, width, wave_length, deltaX, deltaY, fx_max, fy_max):
    # create a mesh_grid to operate in world-coordinates
    x = np.arange(0, height, 1)  # array x
    y = np.arange(0, width, 1)  # array y
    X, Y = np.meshgrid(x - (height / 2), y - (width / 2), indexing='xy')  # meshgrid XY
    fx_0 = height / 2
    fy_0 = width / 2
    k = (2 * math.pi) / wave_length
    theta_x = math.asin((fx_0 - fx_max) * wave_length / (height * deltaX))
    theta_y = math.asin((fy_0 - fy_max) * wave_length / (width * deltaY))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * deltaX) + (math.sin(theta_y) * Y * deltaY)))
    return ref_wave


def numerical_reconstructions(hologram_filter, ref_wave):
    reconstruction = hologram_filter * ref_wave
    phase = np.angle(reconstruction)
    return phase


def autorecons(hologram, wave_length, deltaX, deltaY, diameter):
    hologram, height, width = cut_hologram(hologram)
    fx_max, fy_max, hologram_filter = spatial_filter(hologram, height, width, diameter)
    # ref_wave = reference_wave(height, width, wave_length, deltaX, deltaY, fx_max, fy_max)
    # phase = numerical_reconstructions(hologram_filter, ref_wave)
    peaks = [fx_max, fy_max]
    res = minimize(cost_function, peaks, args=(height, width, hologram_filter, wave_length, deltaX, deltaY),
                   method='TNC', bounds=((fx_max - 2, fx_max + 2), (fy_max - 2, fy_max + 2)), tol=1e-2)
    print(res)
    fx_max = res.x[0]
    fy_max = res.x[1]
    ref_wave = reference_wave(height, width, wave_length, deltaX, deltaY, fx_max, fy_max)
    phase = numerical_reconstructions(hologram_filter, ref_wave)

    return phase

