This cython module computes bessel functions of the first kind (J) for integers order n (0,1,2).
Once it has been installed, you can import it as a standard python module.
Bessel functions of the first kind, denoted as Ja(x), are solutions of Bessel's differential 
equation that are finite at the origin (x = 0) for integer a, and diverge as x approaches zero 
for negative non-integer a. The solution type (e.g.,integer or non-integer) and normalization of 
Ja(x) are defined by its properties below. It is possible to define the function by its Taylor 
series expansion around x = 0.

Bessel functions for python are only supported by SciPy package.
Windows users: please download and install Scipy package.

INSTALLATION REQUIREMENTS

Cython (http://cython.org/#download)

INSTALLATION HOW TO:

pip install Cython
pip install jBessel

or download source and run:

python setup.py install 

Then simply start a Python session and do:
from Bessel import jBessel