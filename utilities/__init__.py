# Importing Image and ImageOps module from PIL package  
from PIL import Image, ImageOps 
#from skimage.io import imsave
import scipy.misc


def salutation():
    print ("Hello world! This is pyDHM library!")
    return

def imageRead (namefile):
        
    Im = Image.open(namefile)
    loadImage = ImageOps.grayscale(Im)
	
    return loadImage

def imageSave (namefile, obj):

    #im = Image.fromarray(obj).convert('RGB')  
    #im.save(namefile)
    #imsave(namefile,obj)
    scipy.misc.imsave(namefile, obj)
	
    return


