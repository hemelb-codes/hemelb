#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# When installed by conda or the default CMake, VTK, ITK, and VMTK do
# not show up to setuptools as they lack dist-info directories. Create
# these.

SITE_PACKAGES_DIR=$(python -c 'import sysconfig; print(sysconfig.get_path("platlib"))')

function make_dist_info() {
    NAME=$1
    if pip list --format json | jq -e ".[] | select (.name == \"$NAME\" )" > /dev/null; then
	echo "pip already knows about $NAME"
	return 0
    fi

    VERSION=$(conda list --json $NAME | jq -r -e '.[].version')
    if [ $? -ne 0 ]; then
	echo "Couldn't get version for $NAME"
	exit 1
    fi
    if [ ! -d $SITE_PACKAGES_DIR/$NAME ]; then
	echo "No package dir '$SITE_PACKAGES_DIR/$NAME'"
	exit 1
    fi
    distinfo=$SITE_PACKAGES_DIR/$NAME-$VERSION.dist-info
    mkdir -p $distinfo
    cat > $distinfo/METADATA <<EOF
Metadata-Version: 2.1
Name: $NAME
Version: $VERSION
EOF
}

make_dist_info vtk
make_dist_info vmtk
