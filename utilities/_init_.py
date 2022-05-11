# Importing Image and ImageOps module from PIL package
from PIL import Image, ImageOps
from matplotlib import pyplot as plt


# Salutation function of the library
def salutation():
    print("Hello world! This is pyDHM library version 1.0")
    return


# Function to read an image file from the disk
def imageRead(namefile):
    Im = Image.open(namefile)
    loadImage = ImageOps.grayscale(Im)

    return loadImage


# Function to display an image
# Inputs:
# inp - The input complex field
# title - The title of the displayed image
def imageShow(inp, title):
    plt.imshow(inp, cmap='gray'), plt.title(title)  # image in gray scale
    plt.show()  # show image

    return
