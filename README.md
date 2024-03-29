pyDHM
=============

An open-source Python library to numerically recover the complex wavefield information of samples from Digital Holographic Microscopy (DHM) for a wide variety of experimental recording configurations. Among others, the library includes:

- Phase retrieval from image-plane and non-image-plane recordings.
- Numerical propagation of wavefields via three different propagators (Angular spectrum approach, Fresnel transform, and Fresnel-Bluestein transform).
- Phase compensation of the tilting angle in off-axis DHM (diffraction-limited and non-diffraction-limited recordings) via three different approaches.
- Phase-shifting methods for both in-line (3-step, 4-step, and 5-step), slightly off-axis acquisitions (quadrature phase-shifting), and two blind phase-shifting methods (No known phase-shifts).
- Phase recovery for systems working in telecentric and non-telecentric (compensation of the quadrature phase factor) regimes.

## Installation

You can install pyDHM from [PyPI](https://pypi.org/project/pyDHM/):

    pip install pyDHM

pyDHM is supported on Python 3.7 and above.

To use pyDHM, you must install opencv (cv2) and scipy in your enviroment.

## Documentation

The pyDHM documentation page can be found in [Project Documentation](https://catrujilla.github.io/pyDHM/).

## How to use

Three sample scripts are provided:

- The “example_numericalPropagation.py” script includes examples of the use of the three numerical propagators implemented in the library to digitally focus complex wavefields.
- The “examples_compensation.py” script includes examples of the use of the three compensation methods for telecentric off-axis DHM holograms and one example of the non-telecentric wavefield compensation.
- The “examples_phaseShitting.py” script includes examples of the use of two phase-shifting retrieval methods for slightly off-axis DHM holograms and three phase-shifting retrieval methods for in-line setups.

More on the use of the pyDHM library can be found in the following academic papers:

Castañeda R, Trujillo C, Doblas A (2022) pyDHM: A Python library for applications in digital holographic microscopy. PLoS ONE 17(10): e0275818.
https://doi.org/10.1371/journal.pone.0275818

R. Castaneda, C. Trujillo, and A. Doblas, "An Open-Source Python library for Digital Holographic Microscopy Imaging," in Imaging and Applied Optics Congress 2022 (3D, AOA, COSI, ISA, pcAOP), Technical Digest Series (Optica Publishing Group, 2022), paper JTh2A.1.
https://opg.optica.org/abstract.cfm?URI=3D-2022-JTh2A.1

## Acknowledgment

The authors would like to thank the 3D Imaging & Display Laboratory, co-lead by Drs. M. Martinez-Corral and G. Saavedra, and Opto-Digital Processing Group, lead by Dr. J. Garcia-Sucerquia for providing the experimental holograms used to validate the pyDHM library.

## About us

The library has been envisioned, designed, and implemented by the Optical Imaging Research laboratory (OIRL) from the University of Memphis and the Applied Optics research group from Universidad EAFIT. The main contributors to the library are:

Dr. Raúl Castañeda-Quintero (University of Memphis) rcstdqnt@memphis.edu

Professor Dr. Ana Doblas (University of Memphis) adoblas@memphis.edu

Professor Dr. Carlos Trujillo (Universidad EAFIT) catrujilla@eafit.edu.co
