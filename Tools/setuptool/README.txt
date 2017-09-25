Step 0 - Requirements
=====================

To build the setuptool you require the following python packages
- numpy
- VTK
- hemeTools (see $HEMELB_ROOT/Tools)
- VMTK (see http://www.vmtk.org/)

And the following non-python software:
- C++11 compiler.
- Boost
- HDF5
- CGAL

To enable the GUI you have to have wx installed and have this working
with your VTK module.

Step 1 is mostly automatic if you have CMake >= 3.2


Step 1 - Configure
==================

You must generate a setup.cfg for the build. Try to use CMake if
possible for this.

Run Cmake for an in-source build. This is usually bad practice, but
required because Python...

> cmake .

Alternatively, copy setup.cfg.in to setup.cfg and edit it appropriately.

Step 2 - Build
==============

Usual python setup.py install

If you're developing, consider instead

> python setup.py build_ext --inplace

Step 3 - Clean
==============

Run clean.sh



