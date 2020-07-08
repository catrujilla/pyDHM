import cv2
from matplotlib import pyplot as plt
import AutoReconstruction

def read_sample():
    # name = 'Basler_acA1920-25um__23117430__20200521_135622465_1.tiff'
    # name = 'Borojo_40x.jpg'
    name = 'DrosophilaMelanogaster.tif'
    # name = 'holotele8x.tif'
    # name = 'holotele10x.tif'
    # name = 'holoRBCtele50x.tif'
    # name = 'image_23.tif'
    hologram = cv2.imread(name, cv2.IMREAD_GRAYSCALE)  # convert automatically to gray scale the image read
    return hologram


def main():
    hologram = read_sample()
    # Borojo
    # phase = AutoReconstruction.autorecons(hologram, 0.633, 6.9, 6.9, 100)
    # DrosophilaMelanogaster
    phase = AutoReconstruction.autorecons(hologram, 0.633, 6.9, 6.9, 100)
    plt.imshow(phase, cmap='gray'), plt.title('phase recontruction')  # image in gray scale
    plt.show()  # show hologram

main()