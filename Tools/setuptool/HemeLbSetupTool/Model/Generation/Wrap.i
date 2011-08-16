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
#include "ConfigGenerator.h"
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

%ignore ConfigGenerator::ClassifySite(Site&);
%include ConfigGenerator.h
%include Iolet.h

%ignore Iterator;
%ignore Vec3::operator%;
%ignore Vec3::operator[];

%extend Vec3 {
  // Python specific wrapping of operator[]
  Vec3::value_type __getitem__(int i) {
    // $self is a pointer to Vec3
    return (*$self)[i];
  }
  void __setitem__(int i, Vec3::value_type val) {
    // $self is a pointer to Vec3
    (*$self)[i] = val;
  }
  %pythoncode %{
    def __str__(self):
        return "Vec3(%s, %s, %s)" % (self[0], self[1], self[2])
  %}
};

%feature("shadow") Vec3::__getitem__ %{
def __getitem__(self, i):
    if i < 0 or i > 2:
        raise IndexError("Vec3 index out of range")
    return $action(self, i)
%}

%feature("shadow") Vec3::__setitem__ %{
def __setitem__(self, i, val):
    if i < 0 or i > 2:
        raise IndexError("Vec3 index out of range")
    return $action(self, i, val)
%}

%ignore Vec3::Magnitude;
%include Index.h
%template (DoubleVector) Vec3<double>;
