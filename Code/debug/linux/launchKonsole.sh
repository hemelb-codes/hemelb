#!/bin/bash
DIR=$(dirname "$0")
debuggerCommandFile=$DIR/resume.gdb
binary=$1
debugger="gdb -q -x $debuggerCommandFile $binary"

shift
# Remove the first to params.
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

