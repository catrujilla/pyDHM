# pyDHM

## An Open-Source Python library for Digital Holographic Microscopy Imaging.

pyDHM is an open-source Python library aimed at Digital Holographic Microscopy (DHM) applications. The pyDHM is a user-friendly library written in the robust programming language of Python that provides a set of numerical processing algorithms for reconstructing amplitude and phase images for a broad range of optical DHM configurations. The pyDHM implements phase-shifting approaches for in-line and slightly off-axis systems and enables phase compensation for telecentric and non-telecentric systems. In addition, pyDHM includes three propagation algorithms for numerical focusing complex amplitude distributions in DHM and digital holography (DH) setups. We have validated the library using numerical and experimental holograms.

### Reconstruction algorithms for different optical DHM configurations 

The selection of the proper computational reconstruction algorithm for a DHM configuration is a critical task to avoid incorrect sample measurements.

![DHM setup](/images/fig1.jpg)

#### How to select the correct reconstruction algorithm?

- Operation of the DHM system based on the interference angle.

![Interference_angle](/images/fig2.JPG)

- Telecentric (d = fTL) vs Non-Telecentric DHM systems (d ≠ fTL)

![Telecentric](/images/fig3.JPG)

- In-focus vs out-of-focus holograms

![In-focus](/images/fig4.JPG)



## pyDHM library structure and examples

The pyDHM library consists of four packages. The first utility package includes essential functions such as reading, displaying images, and filtering Fourier spectrums of holograms. The second package is related to reconstructing in-line and slightly off-axis DHM systems using phase-shifting techniques. The third package reconstructs phase images from off-axis telecentric-based DHM holograms. Finally, the last package includes algorithms for propagating complex amplitude distributions using the angular spectrum, the Fresnel and Fresnel-Bluestein transform approaches. This section explains each package in detail, including the functions and required parameters. We also present sample codes and results for illustrating the performance of the pyDHM library. For more information, please read our pyDHM publication.

### Utility package

The first package in the pyDHM library contains functions for reading and displaying images, computing the Fourier transform (FT), and applying filters to reduce speckle noise. Since the library focuses on DHM applications dealing with complex amplitude distributions, one can display any complex wavefield's amplitude, intensity, or phase map. Although these operations can be straightforwardly implemented in Python for experienced users, this package is aimed to provide compact and user-friendly codes. This package is imported by typing the following code lines:

*from pyDHM import utilities*

The information for each package function is presented bellow, including the declaration statement and the parameters needed. Examples of the use of this package are shown in the upcoming figures. 

#### Available functions in the utility package:

##### imageRead

*imageRead(namefile)*

Function to read an image. The parameter *namefile* corresponds to the name of the image to be opened (e.g., the hologram). 

##### imageShow

*imageShow(inp, name)*

Function to display an image. Two parameters are necessary: *inp* is the data to be visualized (e.g., the loaded hologram, the amplitude, intensity or phase distribution), and *name* is a label for the displayed image.  

##### amplitude

*amplitude(output, log)*

Function to compute the amplitude distribution of the output complex wavefield. Two parameters are necessary: *output* is the complex wavefield distribution, and *log* corresponds to a Boolean variable (e.g., True or False) for applying a common logarithm transformation to the amplitude distribution. 

##### intensity

*intensity(output, log)*

Function to compute the intensity distribution of the output complex amplitude wavefield. Two parameters are necessary: *output* is the complex wavefield distribution, and *log* is the Boolean variable to apply a common logarithm transformation. 

##### phase

*phase(output)*

Function to compute the phase distribution of an output complex wavefield distribution. The only required parameter is *output*, the complex wavefield distribution.

##### Fourier Transform

*FT(input)*

Function to compute the 2D Fourier transform of an image. The only required parameter is the image, *input*.

##### Inverse Fourier Transform

*IFT(input)*

Function to compute the 2D inverse Fourier transform of a spectral image. The only required parameter is the spectral image, *input*.

##### Circular filter

*sfc(field, radius, centX, centY, display)*

Function to filter the Fourier Transform of a hologram using a circular mask. The required parameters are: *field* is the hologram, *radius* is the radius of the circular mask in pixels, and (*centX*, *centY*) are the central pixel positions for the circular mask. The boolean parameter *display*: if True, the spatially filtered hologram is displayed.

##### Rectangular filter

*sfr(field, x1, x2, y1, y2, display)*

Function to filter the Fourier Transform of a hologram using a rectangular mask. The required parameters are: *field* the hologram; (*x1*, *y1*) the pixel coordinates of the upper left corner for the rectangular mask; and (*x2*, *y2*) the pixel coordinates of the lower right corner for the rectangular mask. The boolean parameter *display*: if True, the spatially filtered hologram is displayed.

##### Manual rectangular filter

*sfmr(field, display)*

Function to filter the Fourier Transform of a hologram with a manually selected rectangular mask (OpenCV functionality). The required parameters are *field*, the hologram, and the boolean parameter *display*; if True, the spatially filtered hologram is displayed.

##### Hybrid median-mean

*HM2F(inp, kernel)*

Function to apply the median-mean filter to reduce speckle noise. The parameters are: *inp* the reconstructed amplitude or phase image to be applied the filter; and *kernel* corresponds to the maximum kernel size for the median filter.

### Phase-shifting package

The second package in the pyDHM library contains the phase-shifting strategies for reconstructing the complex amplitude distribution in in-line and slightly off-axis systems. The following code line:

*from pyDHM import phaseShifting*

calls this package. The package is composed of six different phase-shifting (PS) approaches. We have implemented the traditional phase-shifting techniques in which the phase shifts are known using 5 (PS5), 4 (PS4) and 3 (PS3) phase-shifted images. We have also implemented the quadrature PS method (SOSR) and two blind PS approaches using 3 (BPS3) and 2 (BPS2) frames for slightly off-axis DHM systems. The two blind PS approaches require a DHM operating in telecentric regime. The different PS strategies implemented in the package, their definition line statement, and respective parameters are presented bellow. 

#### Available functions in the Phase-shifting package:

##### 5 frames phase-shifting (in-line)

*PS5(inp0, inp1, inp2, inp3, inp4)*

Function to reconstruct the phase distribution using 5 phase-shifted holograms with a phase shift of π/2. The required parameters are the five holograms *inp0* to *inp4*.

##### 4 frames phase-shifting (in-line)

*PS4(inp0, inp1, inp2, inp3)*

Function to reconstruct the phase distribution using 4 phase-shifted holograms with a phase shift equal to π/2. The required parameters are the four holograms *inp0* to *inp3*.

##### 3 frames phase-shifting (in-line)

*PS3(inp0, inp1, inp2)*

Function to reconstruct the phase distribution using 3 phase-shifted holograms with a phase shift equal to π/3. The required parameters are the three holograms *inp* to *inp2*.

##### quadrature method (slightly off-axis)

*SOSR(inp0, inp1, inp2, inp3, upper, wavelength, dx, dy, s=1, steps=4)*

Function to reconstruct the phase distribution using 4 phase-shifted holograms with a phase shift equal to π/2. This quadrature method is based on the SOSR approach proposed by De Nicola et al. [De Nicola S, Ferraro P, Finizio A, Pierattini G. Wave front reconstruction of Fresnel off-axis holograms with compensation of aberrations by means of phase-shifting digital holography. Opt Lasers Eng. 2002;37(4):331–40]. The input parameters are the four holograms, *inp0* to *inp3*; *upper* corresponds to a region for searching the diffraction order; *wavelength* is the illumination wavelength, (*dx*,*dy*) are the pixel sizes along the x- and y- axis, respectively. *s* is a parameter to determine the ROI size in pixels in each dimension to search for the best spatial frequency (size in each dimension equals 1+2s). *steps* is the number of steps for the search inside the ROI in each dimension.

##### Blind 3 raw frames phase-shifting (slightly off-axis)

*BPS3(inp0, inp1, inp2, wavelength, dx, dy)*

Function to reconstruct the phase distribution using 3 phase-shifted holograms with an arbitrary and unknown phase shift. The input parameters are the three holograms (*inp0*-*inp2*); the illumination wavelength (*wavelength*), and the pixel sizes for x and y axes (*dx*, *dy*). This method is valid for slightly off-axis DHM systems operating in telecentric regime. 

##### Blind 2 raw frames phase-shifting (slightly off-axis)

*BPS2(inp0, inp1, wavelength, dx, dy)* 

Function to reconstruct the phase distribution using 2 phase-shifted holograms with an arbitrary and unknown phase shift. This method is valid for slightly off-axis DHM systems operating in telecentric regime. The spectral +1 and -1 components in the spectrum of the recorded hologram should not overlay. The input parameters are the two holograms (*inp0*-*inp1*); the wavelength of the illumination (*wavelength*), and the pixel sizes for x and y axes (*dx*, *dy*).

##### Example of use of the package

![Example 1](/images/fig5.JPG)


### Fully-compensated phase reconstruction package

The third package of the pyDHM library is devoted to the phase reconstruction of DHM holograms without or with minimal perturbations (e.g., fully-compensated reconstructed phase images without distorting sawtooth fringes) using an off-axis system. To call this package you must include the following line:

*from pyDHM import phaseCompensation*

In this package, we have implemented four functions, three functions for holograms recorded in telecentric regime: the full-ROI-search (FRS) function, the efficient-ROI -search (ERS) function, and the cost-function-search (CFS) function. And one function for holograms recorded in non-telecentric regime, the compensation non-telecentric (CNT) function. 

The definition statement and a brief description of each package function is presented bellow.

#### Available functions in the Fully-compensated phase reconstruction package:

##### Full ROI search

*FRS(inp, upper, wavelength, dx, dy, s, step)*

Function to reconstruct fully-compensated phase images using the full ROI search strategy [Trujillo C, Castañeda R, Piedrahita-Quintero P, Garcia-Sucerquia J. Automatic full compensation of quantitative phase imaging in off-axis digital holographic microscopy. Appl Opt. 2016;55(36):10299–306]. The function requires seven parameters: *inp* the hologram, *upper* a boolean variable to determine if the +1 difraction order is located in the upper half of the hologram spectrum, *wavelength* the illumination wavelength, *dx* and *dy* the pixel sizes, and *s* and *step* are parameters for selecting the search region to find the best phase reconstructed image. The default values of these parameters are s = 2, and step = 10. The DHM system should operate in off-axis configuration and telecentric regime.

##### Efficient ROI search

*ERS(inp, upper, wavelength, dx, dy, s, step)*

Function to reconstruct fully-compensated phase images using the efficient ROI search strategy [Obando S, Trujillo C. Computationally efficient phase aberration compensation method for digital holographic microscopy of biological samples. In: Biophotonics Congress 2021. Optical Society of America; 2021]. The function requires seven parameters: *inp* the hologram, *upper* a boolean variable to determine if the +1 difraction order is located in the upper half of the hologram spectrum, *wavelength* the illumination wavelength, *dx* and *dy* the pixel sizes, and *s* and *step* are parameters for selecting the search region to find the best phase reconstructed image. The default values of these parameters are s = 5, and step = 0.2. The DHM system should operate in off-axis configuration and telecentric regime.

##### Cost-function search

*CFS(inp, wavelength, dx, dy)*

Function to reconstruct fully-compensated phase images by minimizing a cost function [Castaneda R, Doblas A. Fast and automatic algorithm to universal recovery of the quantitative phase distribution in digital holographic microscopy. Appl Opt. 2021;60(32):10214–20]. The function requires four parameters: *inp* the hologram, *wavelength* the illumination wavelength, *dx* and *dy* the pixel sizes. The DHM system should operate in off-axis configuration and telecentric regime.

##### Compensation no-telecentric regime

*CNT(inp, wavelength, dx, dy, x1, x2, y1, y2, spatialFilter)*

Function to reconstruct fully-compensated phase images of holograms recorded in non-telecentric regime. The function requires nine parameters: *inp* the hologram, *wavelength* the illumination wavelength, *dx* and *dy* the pixel sizes, and the optional *x1*, *x2*, *y1*, *y2* to input the oocrdinates of a rectangular mask to filter the non-telecentric +1 diffraction order. The units of wavelength, dx, dy should be microns. The DHM system should operate in off-axis configuration and non-telecentric regime. The parameter *spatialFilter* allows to decide whether a manual filter (*spatialFilter == 'sfmr'*) or a rectangular mask filter (*spatialFilter == 'sfr'*, using the x1, x2, y1, y2 parameters) is implemented. This function is based on the proposed method by Kemper et al. [Min J, Yfao B, Ketelhut S, Engwer C, Greve B, Kemper B. Simple and fast spectral domain algorithm for quantitative phase imaging of living cells with digital holographic microscopy. Opt Lett. 2017;42(2):227–30]. We have implemented a search algorithm with nested for loops [Trujillo C, Castañeda R, Piedrahita-Quintero P, Garcia-Sucerquia J. Automatic full compensation of quantitative phase imaging in off-axis digital holographic microscopy. Appl Opt. 2016;55(36):10299–306] to find automatically the best reconstructed phase image (e.g. fully-compensate).

##### Example of use of the package

![Example 2](/images/fig6.JPG)

### Numerical propagation package 

The final package in pyDHM is the numerical propagation package. This package contains the numerical propagation algorithms to compute the scalar complex diffractive wavefield at different propagation distances. The package is called by the following code line:

*from pyDHM import numericalPropagation*

We have implemented three different propagators: angular spectrum *angularSpectrum* the Fresnel transform *fresnel*, and the Fresnel-Bluestein transform *bluestein*. The declaration statement and the parameters needed for each propagator are presented bellow. 

#### Available functions in the Numerical propagation package:

##### Angular spectrum

*angularSpectrum(field, z, wavelength, dx, dy)*

Function to propagate a complex distribution using the angular spectrum approach. *field* is the input complex wavefield to be propagated. The distance to propagate the input wavefield is represented by *z*. *wavelength* is the wavelength of the illumination source used to record the hologram. Finally, *dx* and *dy* are the pixel size for the input and output planes along the x- and y- directions.

##### Fresnel transform 

*fresnel(field, z, wavelength, dx, dy)*

Function to propagate a complex distribution using the Fresnel Transform approach. *field* is the input complex wavefield to be propagated. The distance to propagate the input wavefield is represented by *z*. *wavelength* is the wavelength of the illumination source used to record the hologram. Finally, *dx* and *dy* are the pixel size for the input and output planes along the x- and y- directions.

##### Fresnel-Bluestein transform 

*bluestein(field, z, wavelength, dx, dy, dxout, dyout)*

Function to propagate a complex distribution using the Fresnel-Bluestein transform approach. *field* is the input complex wavefield to be propagated. The distance to propagate the input wavefield is represented by *z*. *wavelength* is the wavelength of the illumination source used to record the hologram. Finally, *dx* and *dy* are the pixel size for the input and output planes along the x- and y- directions. The bluestein propagator function has two additional parameters (e.g., *dxout* and *dyout*) related to the pixel size at the output plane. Note that the pixel size of the input (*dx*, *dy*) and the output (*dxout*, *dyout*) planes can be different in this method.

##### Example of use of the package

![Example 3](/images/fig7.JPG)


### Troubleshooting 

The pyDHM library also incorporates some appropriate error messages to guide the users upon error. In particular, we have created error messages for five particular situations. The first two situations are related to the running environment and having the correct libraries installed. In particular, error messages are displayed if the OpenCV (cv2) and scipy libraries are not installed. On the other hand, an error message is displayed when users incorrectly use one of the functions based on the given holograms. In other words, we have implemented a local function, named *regime*, within the in-line PS techniques (e.g., PS3, PS4 and PS5) and the telecentric-based phase compensation methods (e.g., FRS,ERS, and CFS) that checks the optical configuration of the DHM system (e.g., on/slight off axis versus off axis) based on the spectrum of the input hologram and determines if the called function is suitable for the given hologram. For example, let us assume that the three terms in the hologram's spectrum overlap (e.g., non off-axis configuration) and the user calls the ERS function. Since that function is only valid for holograms recorded in off-axis configuration, an error message (e.g., ‘ERS requires an off-axis hologram’) is displayed. We have also implemented error messages when incorrect values of s and step parameters are used for the SORS, FRS, and ERS functions. Finally, the CNT function enables the spatial filtering of the object frequencies from the hologram spectrum using a popup window (*sfmr*) or a rectangular filter mask defined by the x1, x2, y1, y2 parameters. An error message will appear if the user tries to call the CNT function using the x1, x2, y1, y2 parameters with the option *sfmr*. Also, another error message will appear if the user calls the CNT function using the *sfr* option without inserting the x1, x2, y1, y2 parameters.

In summary, the pyDHM library includes the below errors messages for the following issues:

1) When (x1,y1,x2,y2) values are introduced but the spatial filter option is ‘sfmr’ 

*Please use the sfr option for the spatial filter if (x1,y1,x2,y2) are used as inputs or remove the (x1,y1,x2,y2) as inputs of the CNT function.*

2) When (x1,y1,x2,y2) values are not introduced but the spatial filter option is ‘sfr’ 

*Please use the sfmr option for the spatial filter or insert the values for  (x1,y1,x2,y2) is the spatial filter option is ‘sfr’*

3) If the scipy library is not installed in your desktop/laptop 

*Please install scipy library and import the function minimize*

4) If the OpenCV library is not installed in your desktop/laptop 

*Please install OpenCV library before using the sfmr function and option of the CNT function* 

5) If the hologram used as input of the in-line phase-shifting algorithms (e.g., PS3, PS4, PS5) is recorded in off-axis configuration (e.g., the three terms composing the hologram spectrum do not overlap) 

*PS3/PS4/PS5 algorithm requires on-axis holograms* 

6) If the hologram used as input of the phase compensation algorithms for off-axis telecentric-based DHM systems (e.g., FRS,ERS, CFS) is recorded in on-axis configuration (e.g., the three terms composing the hologram spectrum overlap) 

*FRS/ERS/CFS algorith requires on-axis holograms* 

### Videos: How to use it

#### How to use the angular spectrum function to focus complex distributions computationally

<p align="center">
<iframe width="1519" height="585" src="https://www.youtube.com/embed/3p6Bsh048Hw" title="pyDHM - how to use the angular spectrum function to focus complex distributions computationally" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</p>  

#### How to use the FRS function to reconstruct phase images with minimum phase distortions

<p align="center">
<iframe width="1519" height="585" src="https://www.youtube.com/embed/CMHbF0uoWDk" title="pyDHM - how to use the FRS function to reconstruct phase images with minimum phase distortions" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</p>  

#### How to use the blind 2-frame phase shifting function (BPS2)

<p align="center">
<iframe width="1519" height="585" src="https://www.youtube.com/embed/Z9o0ODe1lUQ" title="pyDHM library - how to use the blind 2-frame phase shifting function (BPS2)" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</p>  

#### How to install pyDHM using pycharm

<p align="center">
<iframe width="1519" height="585" src="https://www.youtube.com/embed/h76nZM6JpXo" title="pyDHM library - how to install using pycharm" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</p>  

### Installation

You can install pyDHM from [PyPI](https://pypi.org/project/pyDHM/):

*pip install pyDHM*

pyDHM is supported on Python 3.7 and above.

To use pyDHM, you must install opencv (cv2) and scipy in your enviroment.


### Funding

This research was partially funded by the Vicerrectoría de Ciencia, Tecnología e Innovación from Universidad EAFIT, and National Science Foundation (NSF) grant number 2042563.


### Credits

pyDHM was developed in python 3.7, numpy 1.18, opencv 4.5 and scipy 1.4.

### Acknowledgment

The authors would like to thank the 3D Imaging & Display Laboratory, co-lead by Drs. M. Martinez-Corral and G. Saavedra, and Opto-Digital Processing Group, lead by Dr. J. Garcia-Sucerquia for providing the experimental holograms used to validate the pyDHM library.

### Citation

If using pyDHM for publication, please kindly cite the following: 

*Castañeda R, Trujillo C, Doblas A (2022) pyDHM: A Python library for applications in digital holographic microscopy. PLoS ONE 17(10): e0275818.*
https://doi.org/10.1371/journal.pone.0275818

*R. Castaneda, C. Trujillo, and A. Doblas, "An Open-Source Python library for Digital Holographic Microscopy Imaging," in Imaging and Applied Optics Congress 2022 (3D, AOA, COSI, ISA, pcAOP), Technical Digest Series (Optica Publishing Group, 2022), paper JTh2A.1.*
https://opg.optica.org/abstract.cfm?URI=3D-2022-JTh2A.1

### Support or Contact 

| Researcher  | email | Google Scholar | 
| ------------- | ------------- |-------------| 
| Raul Castañeda | *sobandov@eafit.edu.co* |  [RaulGoogle](https://scholar.google.com/citations?user=RBtkL1oAAAAJ&hl=es) |
| Ana Doblas| *adoblas@memphis.edu* | [AnaGoogle](https://scholar.google.es/citations?user=PvvDEMYAAAAJ&hl=en) |
| Carlos Trujillo| *catrujilla@eafit.edu.co* | [CarlosGoogle](https://scholar.google.com/citations?user=BKVrl2gAAAAJ&hl=es) |

The Principal Investigators of pyDHM project are Dr. Carlos Trujillo and Dr. Ana Doblas.
