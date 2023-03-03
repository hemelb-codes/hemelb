<!-- This file is part of HemeLB and is Copyright (C) -->
<!-- the HemeLB team and/or their institutions, as detailed in the -->
<!-- file AUTHORS. This software is provided under the terms of the -->
<!-- license in the file LICENSE. -->

# HemeLB geometry generation tool

This tool (and optional GUI) generate HemeLB input geometry files.
This was previously called the "setuptool" so you may see references
to that within the repository.


## Install

Build dependencies:

- Python >= 3.6
- Setuptools >= 21.1
- C++ compiler (C++20)
- CMake >= 3.13
- Pybind11
- scikit-build
- Boost (header only)
- VTK >= 9
- CGAL

Runtime dependencies:
- Python
- Numpy
- PyYAML
- VTK
- VMTK
- wxPython (only for the GUI)

Test dependencies:
- pytest
- hlb (NB: this must be compiled with the same version of numpy as
  used above)


## Conda install

VMTK is hard to build yourself but is easily installed via Conda, so
we present this process here.

### Virtual enviroment setup

We *strongly* recommend that you install into an isolated virtual
environment.

Create the environment for the tool and specify required
packages. They live in the `conda-forge` "channel". The
dependencies are recorded in the `conda-environment.yml` file.

```
conda env create --file conda-environment.yml
```

By default this will use `gmy-tool` as the name, but you can use
anything you wish by adding `--name PREFERED_NAME` to the command
above.


You then need to activate this for your shell session:

```
conda activate gmy-tool
```

### Install

With the environment above, install should be as simple as

```
pip install '.[gui]'
```

(If your pip version is less than 21.3, add the extra flag `--use-feature=in-tree-build`)

If you don't want the GUI, you can drop the `[gui]` extra
specification.

**Note for macOS**: On macOS GUI applications have to be linked
against some special framework. If you don't fix this up you will get
a message:

```
This program needs access to the screen. Please run with a
Framework build of python, and only when you are logged in
on the main display of your Mac.
```

To fix this, you have edit the shebang (`#!`) line in the launcher
script (`hlb-gmy-gui`) to point to an appropriately linked python
executable. This is typically `pythonw` but for conda it is part of
the `python.app` package. Delightfully, this is not in the typical
`$CONDA_PREFIX/bin` directory, instead its location is (2021)
`$CONDA_PREFIX/python.app/Contents/MacOS/python`.

We include a script to fix this for you until pip and scikit build
support doing this automatically:

```
python macos-fix-gui-launcher.py
```

## Install with custom VMTK

The hemelb-codes organisation on GitHub includes a project to build
VMTK  for Ubuntu: https://github.com/hemelb-codes/vmtk-build/

That will hopefully give you an idea for how to proceed. If you are
lucky you can just download the tarball. (Some Ubuntu packages will
have to be installed - see that repo)

## Test

Install pytest the usual way via `pip install pytest`.

The Conda VMTK package forces you to use an old version of numpy. To
ensure that the `hlb` package in `python-tools` is built with the same
one, you need to change to that directory, install Cython (`pip
install cython`) and then install the package with:

```
pip install --no-build-isolation .
```

Run the tests by invoking `py.test` in this directory.


## Run GUI

Ensure your environment is activated then run `hlb-gmy-gui`. There is
basic command line help available:

```
$ hlb-gmy-gui --help
usage: hlb-gmy-gui [-h] [--profile PATH] [--stl PATH] [--geometry PATH]
                   [--xml PATH]

Process an input STL file intosuitable input for HemeLB.

optional arguments:
  -h, --help       show this help message and exit
  --profile PATH   Load the profile to use from an existing file. Other
                   options givenoverride those inthe profile file.
  --stl PATH       The STL file to use as input
  --geometry PATH  Config output file
  --xml PATH       XML output file
```

The terminal will produced a few errors that can ignore, like:
`vtkSTLReader (0x7fdaa773bf10): A FileName must be specified.`. (This
is just VTK trying to display the mesh before the source file is
specified.)


## Profile (.pr2) files

The geometry tool can store the the data it will use to generate a
geometry file. Saving this is highly recommended for reproducibility!

It's a YAML file which can be edited manually. Floating point values
are stored by default in hexadecimal to avoid precision loss
(https://docs.python.org/3/library/stdtypes.html#float.hex) but can be
set in decimal if more convenient.

Paths are relative to the geometry file's location.
