#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
# Radius = 0.75mm, length = 62mm
env PYTHONPATH=../../../../Tools:$PYTHONPATH python -m hemeTools.surfacegenerator.CylinderGenerator 0.75 62 32 poiseuille_flow_test.stl
