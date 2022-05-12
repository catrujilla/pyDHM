import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()


setuptools.setup(
     name='pyDHM',  
     version='1.0',
     scripts=['pyDHM'] ,
     author="Carlos Trujillo",
     author_email="catrujilla@eafit.edu.co",
     description="An open-source Python library to numerically recover the complex wavefield information of samples from Digital Holographic Microscopy (DHM)",
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/catrujilla/pyDHM",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )