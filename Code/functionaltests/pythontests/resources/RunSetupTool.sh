#!/bin/bash
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
env PYTHONPATH=../../../../Tools/setuptool/:$PYTHONPATH ../../../../Tools/setuptool/scripts/hemelb-setup-nogui poiseuille_flow_test.pro
cp poiseuille_flow_test_master.xml poiseuille_flow_test.xml
