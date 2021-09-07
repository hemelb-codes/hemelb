#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vtkPolyData.h>

#include "PybindVTKTypeCaster.h"
//#include <vtkPythonUtil.h>

#include "CylinderGenerator.h"
#include "GeometryGenerator.h"
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
      .def_readwrite("Centre", &Iolet::Centre)
      .def_readwrite("Normal", &Iolet::Normal)
      .def_readwrite("Radius", &Iolet::Radius)
      .def_readwrite("Id", &Iolet::Id)
      .def_readwrite("IsInlet", &Iolet::IsInlet);

  // ABC GeometryGenerator
  auto gmy_gen =
      py::class_<GeometryGenerator>(mod, "GeometryGenerator")
          .def("GetOutputGeometryFile",
               &GeometryGenerator::GetOutputGeometryFile)
          .def("SetOutputGeometryFile",
               &GeometryGenerator::SetOutputGeometryFile)
          .def("GetIolets", py::overload_cast<>(&GeometryGenerator::GetIolets))
          .def("SetIolets", &GeometryGenerator::SetIolets)
          .def("SetOriginWorking", &GeometryGenerator::SetOriginWorking)
          .def("SetSiteCounts", &GeometryGenerator::SetSiteCounts)
          .def("Execute", &GeometryGenerator::Execute);

  // Concrete PolyDataGenerator
  py::class_<PolyDataGenerator>(mod, "PolyDataGenerator", gmy_gen)
      .def(py::init<>())
      .def("SetSeedPointWorking", py::overload_cast<double, double, double>(
                                      &PolyDataGenerator::SetSeedPointWorking))
      .def("GetClippedSurface", &PolyDataGenerator::GetClippedSurface)
      .def("SetClippedSurface", &PolyDataGenerator::SetClippedSurface);

  // Cylinder
  py::class_<CylinderGenerator>(mod, "CylinderGenerator", gmy_gen)
      .def(py::init<>())
      .def("SetCylinderCentre", &CylinderGenerator::SetCylinderCentre)
      .def("SetCylinderAxis", &CylinderGenerator::SetCylinderAxis)
      .def("SetCylinderRadius", &CylinderGenerator::SetCylinderRadius)
      .def("SetCylinderLength", &CylinderGenerator::SetCylinderLength);
}
