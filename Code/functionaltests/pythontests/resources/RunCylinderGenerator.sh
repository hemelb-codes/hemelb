#!/bin/bash

# Radius = 0.75mm, length = 62mm
env PYTHONPATH=../../../../Tools:$PYTHONPATH python -m hemeTools.surfacegenerator.CylinderGenerator 0.75 62 32 poiseuille_flow_test.stl
