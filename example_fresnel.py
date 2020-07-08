import cv2
from matplotlib import pyplot as plt
import fresnel


def read_sample():
    namefile = 'mask.jpg'
    #namefile = 'die_1.jpg'
    hologram = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
    plt.imshow(hologram, cmap='gray'), plt.title('hologram')  # image in gray scale
    plt.show()  # show hologram
    return hologram


def main():
    hologram = read_sample()
    # mask
    intensitie = fresnel.fr(hologram, 632.8, 0.05, 5.2, 5.2)
    # die
    #intensitie = fresnel.fr(hologram, 632.8, 1.05, 5.2, 5.2)
    plt.imshow(intensitie, cmap='gray'), plt.title('phase recontruction')  # image in gray scale
    plt.show()  # show hologram


main()