#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
DIR=$(dirname "$0")
debuggerCommandFile=$DIR/resume.gdb
binary=$2
debugger="gdb -q -x $debuggerCommandFile $binary"

shift
# Remove the first to params.
# Positional params now hold the process ids we
# want to attach to.
    
command="gnome-terminal"
rank=0

until [ -z "$1" ] # all params used up
do
    pid=$1
    shift
    
    if [ "$rank" -eq "0" ]; then
	command="${command} --window -e '${debugger} $pid'"
    else
	command="${command} --tab -e '${debugger} $pid'"
    fi
    
    ((rank++))
done

echo $command
eval $command

