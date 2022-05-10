// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vtkPolyData.h>

#include "PybindVTKTypeCaster.h"

#include "Index.h"
#include "Iolet.h"
#include "PolyDataGenerator.h"

namespace {
using DV = hemelb::util::Vector3D<double>;
}

namespace py = pybind11;

PYBIND11_MODULE(Generation, mod) {
  mod.doc() = "Wrapped geometry generation code";

  // Map Vector3D exceptions into Python IndexError
  py::register_exception<IndexError>(mod, "DVIndexError", PyExc_IndexError);

  // Expose the required bits of Vector3D<double>
  py::class_<DV>(mod, "DoubleVector")
      // __init__
      .def(py::init(&DV::Zero))                 // Zero
      .def(py::init<double>())                  // All the same
      .def(py::init<double, double, double>())  // Three components
      // Attribute acces to elements
      .def_readwrite("x", &DV::x)
      .def_readwrite("y", &DV::y)
      .def_readwrite("z", &DV::z)
      // get/set via []
      .def("__getitem__", [](DV const& self, int i) { return self[i]; })
      .def("__setitem__", [](DV& self, int i, double v) { self[i] = v; })
      .def(py::self / double());

  // Now the Iolet
  py::class_<Iolet>(mod, "Iolet")
      .def(py::init<>())
      .def_readwrite("Id", &Iolet::Id)
      .def_readwrite("IsInlet", &Iolet::IsInlet);

  // PolyDataGenerator
  py::class_<PolyDataGenerator>(mod, "PolyDataGenerator")
      .def(py::init<>())
      .def("SetSeedPointWorking", py::overload_cast<double, double, double>(
                                      &PolyDataGenerator::SetSeedPointWorking))
      .def("GetClippedSurface", &PolyDataGenerator::GetClippedSurface)
      .def("SetClippedSurface", &PolyDataGenerator::SetClippedSurface)
      .def("GetOutputGeometryFile", &PolyDataGenerator::GetOutputGeometryFile)
      .def("SetOutputGeometryFile", &PolyDataGenerator::SetOutputGeometryFile)
      .def("SetIolets", &PolyDataGenerator::SetIolets)
      .def("SetOriginWorking", &PolyDataGenerator::SetOriginWorking)
      .def("SetNumberOfLevels", &PolyDataGenerator::SetNumberOfLevels)
      .def("SetTriangleLevel", &PolyDataGenerator::SetTriangleLevel)
      .def("Execute", &PolyDataGenerator::Execute);
  ;
}
