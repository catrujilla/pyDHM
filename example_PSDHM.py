import cv2
from matplotlib import pyplot as plt
import numpy as np
import PSDHBlind


def read_sample():
    # Neuron
    #image1 = 'NeuronSilhoutte_1.tif'
    #image2 = 'NeuronSilhoutte_2.tif'
    #image3 = 'NeuronSilhoutte_3.tif'
    # Fresnel lens
    image1 = '1_19.tiff'
    image2 = '2_19.tiff'
    image3 = '3_19.tiff'
    # UN sample
    #image1 = 'UNOD-Test_1.tif'
    #image2 = 'UNOD-Test_2.tif'
    #image3 = 'UNOD-Test_3.tif'

    image1 = cv2.imread(image1, cv2.IMREAD_GRAYSCALE)
    image2 = cv2.imread(image2, cv2.IMREAD_GRAYSCALE)
    image3 = cv2.imread(image3, cv2.IMREAD_GRAYSCALE)
    return image1, image2, image3


def image_show(image, title):
    plt.imshow(image, cmap='gray'), plt.title(title)  # image in gray scale
    plt.show()  # show hologram

def main():
    image1, image2, image3 = read_sample()
    # neuron
    #phase3 = PSDHBlind.raw3(image1, image2, image3, 0.405, 6.9, 6.9)
    #phase2 = PSDHBlind.raw2(image1, image2, 0.405, 6.9, 6.9)
    # Fresnel lens
    phase3 = PSDHBlind.raw3(image1, image2, image3, 0.633, 6.9, 6.9)
    phase2 = PSDHBlind.raw2(image1, image2, 0.633, 6.9, 6.9)
    # UN sample
    #phase3 = PSDHBlind.raw3(image1, image2, image3, 0.405, 6.9, 6.9)
    #phase2 = PSDHBlind.raw2(image1, image2, 0.405, 6.9, 6.9)
    image_show(phase3, 'phase 3 raw')
    image_show(phase2, 'phase 2 raw')

main()