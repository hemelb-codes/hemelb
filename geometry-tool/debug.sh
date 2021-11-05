#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# Append '-g' and '-O0' to the compile flags and build the extension
# in place (i.e. $ python setup.py build_ext --inplace)
# Then treat this script just like the python executable
# (e.g. ./debug.sh scripts/hemelb-setup-nogui /path/to/profile.pro)

# Create a file break.gdb and add a series of GDB break statements.
# They will be set before execution starts.

cat > prebreak.gdb <<EOF
set breakpoint pending on
EOF

cat > run.gdb <<EOF
run $@
continue
EOF

export PYTHONPATH=$CWD:$PYTHONPATH
if [ -e break.gdb ]; then
    break="-x prebreak.gdb -x break.gdb"
fi

gdb $break -x run.gdb python
