import cv2
		
def salutation():
    print ("Hello world! This is pyDiffraction library!")
    return

def imread (namefile):
        
    loadImage = cv2.imread(namefile, cv2.IMREAD_GRAYSCALE)  
	
    return loadImage
