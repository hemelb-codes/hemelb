#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
DIR=$(dirname "$0")
debuggerCommandFile=$DIR/resume.gdb

binary=$1
shift

# Build the debugger invocation
debugger="gdb -q -x $debuggerCommandFile"
# If the -file flag is given then also run those GDB commands
if [ "$1" == "-file" ]; then
  debugger="$debugger -x $2"
  shift 2
fi
# Add the binary
debugger="$debugger $binary"

# Positional params now hold the process ids we
# want to attach to.

tabsFile="debug.tabs"
command="konsole --tabs-from-file $tabsFile"
rank=0

until [ -z "$1" ] # all params used up
do
    pid=$1
    shift
    
    if [ "$rank" -eq "0" ]; then
	    echo "title: Rank $rank ;; command: $debugger $pid ;; workdir: ~" > $tabsFile
    else
	    echo "title: Rank $rank ;; command: $debugger $pid ;; workdir: ~" >> $tabsFile
	fi
    
    ((rank++))
done

echo $command
eval $command

sleep 1
rm $tabsFile

