#!/bin/bash

debuggerCommandFile=$1
binary=$2
debugger="gdb -q -x $debuggerCommandFile $binary"

shift 2 
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

