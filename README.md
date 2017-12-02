# README #

pyres is a Python library for managing and inverting electrical resistivity data using the open source R2 and R3 (ongoing development in pyres) packages from Andrew Binley's website:

http://www.es.lancs.ac.uk/people/amb/Freeware/freeware.htm

### Python package for modeling electrical resistivity data ###

* pyres is a Python wrapper for R2 and R3 allowing mesh creation, forward and inverse modeling, and simple output visualization
* Version 1.0

### Citation ###

Befus, K.M. (2017), pyres: A Python Wrapper for Electrical Resistivity Modeling with R2, J. Geophys. Eng., doi: 10.1088/1742-2140/aa93ad.

### Installation and testing ###

* Install using setup.py
* Or with pip: "pip install pyresistivity"
* Dependencies: numpy, scipy, and matplotlib beyond default python packages
* How to run tests: Test with R2_examples\Surface_*.py scripts

### Examples of pyres use ###

Provided with the pyres code are numerous examples of how to use pyres. These examples correspond to the examples in the R2 documentation and are explained therein. One example of using R2 with field data is also included. 

The main R2 example files are named Surface_1 through Surface_8 with an additional descriptor at the end of the script's name for the resistivity array type used (e.g., Wenner or dpdp for dipole-dipole surveys).

The default installation directory structure is for Windows:

	C:\ER             # main directory holding all data, R2, gmsh, and example scripts.
		\R2_examples  # directory with example scripts and default output location for scripts
		\R2           # R2 installation directory with default structure
		\gmsh-version # gmsh directory

	
### For more information ###

* Please contact Kevin Befus to report bugs or request assistance with pyres: kbefus@uwyo.edu (however time for providing support is very limited)
* Please refer to the documentation at http://www.es.lancs.ac.uk/people/amb/Freeware/freeware.htm for assistance with R2 and R3.
* Please refer to http://gmsh.info for assistance with gmsh.