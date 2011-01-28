#!/bin/bash

# Call with one argument: the machine to generate a build script for.

runDir=$PWD
makeDir=$PWD/`dirname $0`
rootDir=$makeDir/..

if [ -f $makeDir/Makefile.$1 ]; then
    
    cd $rootDir
    make clean_all > /dev/null
    
    export HEMELB_MACHINE=$1
    make -n | \
	sed 's/\/Users\/rupert\/working\/hemelb\/Code\///g;' | \
	sed 's/Rules generated.../#!\/bin\/bash/;' \
	>$runDir/build_on_$HEMELB_MACHINE.sh
    
    chmod a+x $runDir/build_on_$HEMELB_MACHINE.sh
else
    echo "Could not find Makefile for machine $1 in $makeDir!"
fi