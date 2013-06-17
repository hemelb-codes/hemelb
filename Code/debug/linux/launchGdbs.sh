#!/bin/bash

DIR=$(cd $(dirname "$0"); pwd)
# Try to guess which environment to run appropriate debugger script
if [ "$(pidof ksmserver)" ]; then
	# KDE
	$DIR/launchKonsole.sh "$@"
elif [ "$(pidof gnome-session)" ]; then
    # GNOME
    $DIR/launchGnomeTerminal.sh "$@"
else
	# Panic
	echo "Can't figure out environment"
fi
