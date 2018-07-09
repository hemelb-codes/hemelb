#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

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
