pyDHM
=============

An open-source Python library to numerically recover the complex wavefield information of samples from Digital Holographic Microscopy (DHM) for a wide variety of experimental recording possibilities. Among others, the library includes:

- Phase retrieval from image-plane and non-image-plane recordings.
- Numerical propagation of wavefields via three different propagators (Angular spectrum approach, Fresnel transform, and Fresnel-Bluestein transform).
- Phase compensation of the tilting angle in off-axis DHM (diffraction-limited and non-diffraction-limited recordings) via three different approaches.
- Phase-shifting methods for both in-line (3-step, 4-step, and 5-step), slightly off-axis acquisitions (quadrature phase-shifting), and two blind phase-shifting methods (No known phase-shifts).
- Phase recovery for systems working in telecentric and non-telecentric (compensation of the quadrature phase factor) regimes.

Three sample scripts are provided using the library functions:

- The “example_numericalPropagation.py” script includes examples of the use of the three numerical propagators implemented in the library to focus DHM reconstructions.
- The “examples_compensation.py” script includes examples of the use of the three compensation methods for telecentric off-axis DHM holograms and one example of the non-telecentric wavefield compensation.
- The “examples_phaseShitting.py” script includes examples of the use of two phase-shifting retrieval methods for slightly off-axis DHM holograms and three phase-shifting retrieval methods for in-line setups.

The library has been envisioned, designed, and implemented by the Optical Imaging Research laboratory (OIRL) from the University of Memphis and the Applied Optics research group of Universidad EAFIT. The main contributors to the library are:

PhD candidate Raúl Castañeda-Sepulveda (University of Memphis)

Professor Dr. Ana Doblas (University of Memphis)

Professor Dr. Carlos Trujillo (Universidad EAFIT)
