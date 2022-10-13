import pathlib
from setuptools import find_packages, setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text(encoding="utf8")

# This call to setup() does all the work
setup(
    name="pyDHM",
    version="1.0.5",
    description="An open-source Python library to numerically recover the complex wavefield information of samples from Digital Holographic Microscopy (DHM)",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/catrujilla/pyDHM",
    author="Carlos Trujillo",
    author_email="catrujilla@eafit.edu.co",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["pyDHM"],
    include_package_data=True,
    install_requires=["matplotlib", "scipy"],
)