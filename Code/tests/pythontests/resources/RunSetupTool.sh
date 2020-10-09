#!/bin/bash
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
env PYTHONPATH=../../../../Tools/setuptool/:$PYTHONPATH ../../../../Tools/setuptool/scripts/hemelb-setup-nogui poiseuille_flow_test.pro
cp poiseuille_flow_test_master.xml poiseuille_flow_test.xml
