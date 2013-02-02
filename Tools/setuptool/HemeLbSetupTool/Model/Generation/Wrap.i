// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 
%module Generation

%include "cpointer.i"
%pointer_class(Iolet, IoletPtr)

%include "std_string.i"
%include "std_vector.i"
namespace std {
  //%template(DoubleVector) vector<double>;
  %template(IoletPtrVector) vector<Iolet*>;
}
// This is needed declare Iolet before the std::vector support code. But note that we must include Python.h before anything else
%begin %{
#include <Python.h>
#include "Iolet.h"
#include <stddef.h> 
%}

%{
#include "Index.h"
#include "GeometryGenerator.h"
#include "CylinderGenerator.h"
#include "SquareDuctGenerator.h"
#include "PolyDataGenerator.h"
#include "vtkPolyData.h"
#include "vtkOBBTree.h"
#include <sstream>
%}

// TODO: remove the copy-pasta for the 2 vtk classes
%typemap (in) vtkPolyData* {
  int fail = 0;
  // New reference
  PyObject* thisObj = PyObject_GetAttrString($input, "__this__");
  if (thisObj == NULL){
    fail = 1;
  } else {
    // Convert from a Py str
    char* cstr = PyString_AsString(thisObj);
    if (cstr == NULL) {
      fail = 1;
    } else {
      std::string thisStr(cstr);
      // Get the part of the string between the 1st 2 underscores; this holds the pointer part
      size_t under1 = thisStr.find('_');
      size_t under2 = thisStr.find('_', under1 + 1);
      std::string ptrStr = thisStr.substr(under1+1, under2-under1-1);
      // Convert to an unsigned int of the right type.
      std::stringstream ss;
      ss << std::hex << ptrStr;
      uintptr_t ptrUint;
      ss >> ptrUint;
      // Cast to the right type
      $1 = ($1_type)ptrUint;
      // Cancel the reference we took
      Py_DECREF(thisObj);
    }
  }
  if (fail) {
    Py_XDECREF(thisObj);
    SWIG_exception_fail(SWIG_ArgError(SWIG_ERROR), "in method '" "$symname" "', argument " "$argnum" " of type '" "$1_type""'");
  }
 }%typemap (in) vtkOBBTree* {
  int fail = 0;
  // New reference
  PyObject* thisObj = PyObject_GetAttrString($input, "__this__");
  if (thisObj == NULL){
    fail = 1;
  } else {
    // Convert from a Py str
    char* cstr = PyString_AsString(thisObj);
    if (cstr == NULL) {
      fail = 1;
    } else {
      std::string thisStr(cstr);
      // Get the part of the string between the 1st 2 underscores; this holds the pointer part
      size_t under1 = thisStr.find('_');
      size_t under2 = thisStr.find('_', under1 + 1);
      std::string ptrStr = thisStr.substr(under1+1, under2-under1-1);
      // Convert to an unsigned int of the right type.
      std::stringstream ss;
      ss << std::hex << ptrStr;
      uintptr_t ptrUint;
      ss >> ptrUint;
      // Cast to the right type
      $1 = ($1_type)ptrUint;
      // Cancel the reference we took
      Py_DECREF(thisObj);
    }
  }
  if (fail) {
    Py_XDECREF(thisObj);
    SWIG_exception_fail(SWIG_ArgError(SWIG_ERROR), "in method '" "$symname" "', argument " "$argnum" " of type '" "$1_type""'");
  }
 }
 
%typemap(throws) GenerationError %{
  PyErr_SetString(PyExc_RuntimeError, $1.what());
  SWIG_fail;
%}
%include GeometryGenerator.h
//%ignore PolyDataGenerator::ClassifySite(Site&);
%include PolyDataGenerator.h
//%ignore CylinderGenerator::ClassifySite(Site&);
%include CylinderGenerator.h
%include SquareDuctGenerator.h
%include Iolet.h

// Raise a python IndexError when we get a C++ IndexError
%typemap(throws) IndexError %{
  PyErr_SetNone(PyExc_IndexError);
  SWIG_fail;
%}
// Mark these as catching IndexError
%catches(IndexError) hemelb::util::Vector3D::__getitem__;
%catches(IndexError) hemelb::util::Vector3D::__setitem__;

%ignore IndexError;
%ignore hemelb::util::Vector3DIterator;
%ignore hemelb::util::Vector3DBase;
%ignore hemelb::util::Vector3D::operator%;
%ignore hemelb::util::Vector3D::operator[];

%extend hemelb::util::Vector3D {
  // Python specific wrapping of operator[]
  hemelb::util::Vector3D::value_type __getitem__(int i) {
    // $self is a pointer to Vector3D
    return (*$self)[i];
  }
  void __setitem__(int i, hemelb::util::Vector3D::value_type val) {
    // $self is a pointer to Vector3D
    (*$self)[i] = val;
  }
  %pythoncode %{
    def __str__(self):
        return "Vector3D(%s, %s, %s)" % (self[0], self[1], self[2])
  %}

  // Multiplication by a scalar
  hemelb::util::Vector3D __mul__(hemelb::util::Vector3D::value_type divisor) {
    return (*$self) * divisor;
  }
  hemelb::util::Vector3D& __imul__(hemelb::util::Vector3D::value_type divisor) {
  	(*$self) *= divisor;
    return (*$self);
  }
  
  // Division by a scalar
  hemelb::util::Vector3D __div__(hemelb::util::Vector3D::value_type divisor) {
    return (*$self) / divisor;
  }
  hemelb::util::Vector3D& __idiv__(hemelb::util::Vector3D::value_type divisor) {
  	(*$self) /= divisor;
    return (*$self);
  }
};

%feature("shadow") hemelb::util::Vector3D::__getitem__ %{
def __getitem__(self, i):
    if i < 0 or i > 2:
        raise IndexError("Vector3D index out of range")
    return $action(self, i)
%}

%feature("shadow") hemelb::util::Vector3D::__setitem__ %{
def __setitem__(self, i, val):
    if i < 0 or i > 2:
        raise IndexError("Vector3D index out of range")
    return $action(self, i, val)
%}

%ignore hemelb::util::Vector3D::GetMagnitude;
// Swig doesn't follow #includes, so we must manuall %include Vector3D.h
%include util/Vector3D.h
%include Index.h
%template (DoubleVector) hemelb::util::Vector3D<double>;
