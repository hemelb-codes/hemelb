#!/usr/bin/env python

## Program:   PyNS
## Module:    Bessel.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

cdef extern from "complexobject.h":
    struct Py_complex:
        double real
        double imag
    ctypedef class __builtin__.complex [object PyComplexObject]:
        cdef Py_complex cval

def Bessel(int n, double complex arg):
    '''
    Compute Bessel functions of the first kind (J) for integers order n (0,1,2).
    This function takes as input the integer order n (0,1,2) and the argument of the function as a complex number a+bj.
    Bessel functions of the first kind, denoted as Ja(x), are solutions of Bessel's differential equation that are finite 
    at the origin (x = 0) for integer a, and diverge as x approaches zero for negative non-integer a. The solution type 
    (e.g., integer or non-integer) and normalization of Ja(x) are defined by its properties below. It is possible to define 
    the function by its Taylor series expansion around x = 0
    '''
    
    cdef double complex z, zproduct, zanswer, zarg
    cdef double k
    cdef int i
    if n > 2 or n < 0:
        import sys
        sys.exit('Error, index n has to be 0, 1 or 2')
    z = 1. + 0.j
    zproduct = 1. + 0.j
    zanswer = 1. + 0.j
    zarg = -0.25 * (arg * arg)
    
    for i in range(0, 1000):
        k = (i+1.)*(i+1.+n)
        z = (1./(k))*(z*zarg)
        if abs(z) < 1e-20:
            break
        zanswer = zanswer + z
    for i in range(0,n):
        zproduct = zproduct * 0.5 * arg
    zanswer = zanswer * zproduct
    return zanswer