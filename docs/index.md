# pyDHM

## An Open-Source Python library for Digital Holographic Microscopy Imaging.

pyDHM is an open-source Python library aimed at Digital Holographic Microscopy (DHM) applications. The pyDHM is a user-friendly library written in the robust programming language of Python that provides a set of numerical processing algorithms for reconstructing amplitude and phase images for a broad range of optical DHM configurations. The pyDHM implements phase-shifting approaches for in-line and slightly off-axis systems and enables phase compensation for telecentric and non-telecentric systems. In addition, pyDHM includes three propagation algorithms for numerical focusing complex amplitude distributions in DHM and digital holography (DH) setups. We have validated the library using numerical and experimental holograms.

### Reconstruction algorithms for different optical DHM configurations 

The selection of the proper computational reconstruction algorithm for a DHM configuration is a critical task to avoid incorrect sample measurements.

![DHM setup](/images/fig1.jpg)

How to select the correct reconstruction algorithm?

- Operation of the DHM system based on the interference angle.

![Interference_angle](/images/fig2.JPG)

- Telecentric (d = fTL) vs Non-Telecentric DHM systems (d ≠ fTL)

![Telecentric](/images/fig3.JPG)

- In-focus vs out-of-focus holograms

![In-focus](/images/fig4.JPG)



## pyDHM library structure and examples

The pyDHM library consists of four packages. The first utility package includes essential functions such as reading, displaying images, and filtering Fourier spectrums of holograms. The second package is related to reconstructing in-line and slightly off-axis DHM systems using phase-shifting techniques. The third package reconstructs phase images from off-axis telecentric-based DHM holograms. Finally, the last package includes algorithms for propagating complex amplitude distributions using the angular spectrum, the Fresnel and Fresnel-Bluestein transform approaches. This section explains each package in detail, including the functions and required parameters. We also present sample codes and results for illustrating the performance of the pyDHM library. For more information, please read our pyDHM publication.

### Utility package

The first package in the pyDHM library contains functions for reading and displaying images, computing the Fourier transform (FT), and applying filters to reduce speckle noise (62). Since the library focuses on DHM applications dealing with complex amplitude distributions, one can display any complex wavefield's amplitude, intensity, or phase map. Although these operations can be straightforwardly implemented in Python for experienced users, this package is aimed to provide compact and user-friendly codes. This package is imported by typing the following code lines, from pyDHM import utilities. Table 1 shows the information for each package function, including the declaration statement and the parameters needed. Examples of the use of this package are shown in the upcoming figures. 

#### Available functions in the utility package:

##### imageRead

imageRead(namefile)

Function to read an image. The parameter namefile corresponds to the name of the image to be opened (e.g., the hologram). 

##### imageShow

imageShow(inp, name)

Function to display an image. Two parameters are necessary: inp is the data to be visualized (e.g., the load hologram, the amplitude, intensity or phase distribution), and a name is a label for the displayed image.  


<p align="center">
<iframe width="560" height="315" src="https://www.youtube.com/embed/sUeVBAqYXJU" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</p>  


Click here for download the DRF code for MATLAB. 
* [Download MATLAB script](https://drive.google.com/drive/folders/1su6cW7JX1s3KXNQFdRl8nBjlIo2jobZ-?usp=sharing)


### Sample holograms

Here are the hologramas from where the videos shown prevously were exctracted. You can download the folders with the whole set of holograms and try them yourself.

- [USAF Test](https://drive.google.com/drive/folders/1bAfjlpGmKy6tdngCz6NAa601YNRNC2S5?usp=sharing)
- [Red blood cells](https://drive.google.com/drive/folders/1McD-yl8pHNrWxM6Vdd6p0wdMfufiy3BE?usp=sharing)
- [Spoiled vinaigrette](https://drive.google.com/drive/folders/1j4XBxeqpnIbFAP7NSf7gWLpAC7lVeAIX?usp=sharing)


### Funding
This project has received funding from EAFIT University (Medellin, Antioquia, Colombia).

### Citation
If using this information for publication, please kindly cite the following paper:

Obando-Vásquez, S (2021). "Implementación de un algoritmo para recuperación de información de fase de muestras biológicas dinámicas a través de microscopía holográfica digital"[Unpublished Bachelor's thesis] Universidad EAFIT. 

Carlos Trujillo, Raúl Castañeda, Pablo Piedrahita-Quintero, and Jorge Garcia-Sucerquia, "Automatic full compensation of quantitative phase imaging in off-axis digital holographic microscopy," Appl. Opt. 55, 10299-10306 (2016)

### Support or Contact 

| Researcher  | email | Google Scholar | 
| ------------- | ------------- |-------------| 
| Sofia Obando-Vasquez | *sobandov@eafit.edu.co* |  | 
| Ana Doblas| *adoblas@memphis.edu* | [AnaGoogle](https://scholar.google.es/citations?user=PvvDEMYAAAAJ&hl=en) |
| Carlos Trujillo| *catrujilla@eafit.edu.co* | [CarlosGoogle](https://scholar.google.com/citations?user=BKVrl2gAAAAJ&hl=es) |

This Research is a collaboration between Doblas’ and Trujillos’ research lab.
