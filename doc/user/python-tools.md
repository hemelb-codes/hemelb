# HemeLB python tools

These tools are a set of Python libraries and some command-line tools
that use them, for tasks such as converting HemeLB input and output
files into other formats.

Dependencies:

- Python >= 3.6
- setuptools
- Cython
- Numpy
- VTK (python bindings)
- C compiler

## Install

1. Create and activate a virtual environment
2. Install with pip:
```bash
cd hemelb/python-tools
pip install .
```

## Use

For now, please run the script with the `-h` option for help.

Available command line tools:

- `hlb-gmy-decompress`: decompress a geometry file.

- `hlb-gmy-compress`: compress an uncompressed geometry file.

- `hlb-dump-extracted-properties`: convert an extraction file to CSV.

- `hlb-gmy-selfconsistent`: check if a geometry file is self-consistent.

- `hlb-gmy-countsites`: print some basic information about how many sites in a geometry file


Runnable modules (run with `python -m`):

Convert an XML + geometry file to VTK unstructured grid (.vtu)
```
python -m hlb.converters.GmyUnstructuredGridReader path/to/config.xml
```

Convert an extracted property file to VTK unstructured grid (.vtu).
Note you also require EITHER the .xml/.gmy used for the run OR the output of `GmyUnstructuredGridReader` to provide the
geometry information needed:
```
python -m hlb.converters.ExtractedPropertyUnstructuredGridReader geometry.vtu data.xtr
```

