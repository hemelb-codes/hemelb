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
- C++ compiler
- CMake
- Pybind11
- scikit-build
- Boost (header only)
- VTK
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
- hemeTools (NB: this must be compiled with the same version of
  numpy as used above)


## Conda install

VMTK is hard to build yourself but is easily installed via Conda, so
we present this process here.

VMTK requires VTK 8 (unfortunately it does not yet work with VTK 9,
required by the main HemeLB application...).

### Virtual enviroment setup

We *strongly* recommend that you install into an isolated virtual
environment.

Create the environment for the tool and specify required packages. The
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

The VTK and VMTK packages installed by conda do not appear to
setuptools as they don't have the proper metadata included. You can
run the script `bodge-packages-for-setuptools.sh` to fix this.

### Install

With the environment above, install should be as simple as

```
pip install --use-feature=in-tree-build '.[gui]'
```

At some point setuptools will not require the extra flag, but for now
it does.

If you don't want the GUI, you can drop the `[gui]` extra
specification.

### Test

Install pytest the usual way via `pip install pytest`.

The Conda VMTK package forces you to use an old version of numpy. To
ensure that the hemeTools package is built with the same one, you need
to change to that directory, install Cython (`pip install cython`) and
then install the package with:

```
pip install --no-build-isolation .
```

Run the tests by invoking `py.test` in this directory.
