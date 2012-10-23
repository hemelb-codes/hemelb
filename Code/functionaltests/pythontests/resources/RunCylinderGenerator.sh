#!/bin/bash
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
# Radius = 0.75mm, length = 62mm
env PYTHONPATH=../../../../Tools:$PYTHONPATH python -m hemeTools.surfacegenerator.CylinderGenerator 0.75 62 32 poiseuille_flow_test.stl
