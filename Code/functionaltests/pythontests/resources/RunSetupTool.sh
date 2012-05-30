#!/bin/bash

env PYTHONPATH=../../../../Tools/setuptool/:$PYTHONPATH ../../../../Tools/setuptool/scripts/hemelb-setup-nogui poiseuille_flow_test.pro
cp poiseuille_flow_test_master.xml poiseuille_flow_test.xml
