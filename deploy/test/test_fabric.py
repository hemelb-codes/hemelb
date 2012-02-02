#!/usr/bin/env python
# encoding: utf-8
"""
test_machine_environment.py

Created by James Hetherington on 2012-01-19.
Copyright (c) 2012 UCL. All rights reserved.
"""
import unittest
import sys
import copy

from ..fab import *

class TestFabric(unittest.TestCase):
    def setUp(self):
    	#Update the user config with testing example
    	user_config=yaml.load(open(os.path.join(env.localroot,'deploy','test','machines_user.yml')))
    	execute(planck) #Default machine target is assumed as planck.
    	#Monkeypatch the fabric commands to do nothing, but record what they would have done
    	sys.modules['deploy.fab'].run=lambda command: self.commands.append(command)
    	sys.modules['deploy.fab'].local=lambda command: self.commands.append(command)
    	sys.modules['deploy.fab'].put=lambda source,target: self.commands.append("put "+source+" "+target)
    	sys.modules['deploy.fab'].rsync_project=lambda **args: self.commands.append("rsync "+args['local_dir']+" "+args['remote_dir'])
    	sys.modules['deploy.fab'].generate=lambda profile: self.commands.append("generate %g %g %g"%(profile.VoxelSize, profile.Steps , profile.Cycles) )
    	self.commands=[]
    	env.build_number='abcd1234'
    def assertCommandCount(self,should_be):
        self.assertEqual(len(self.commands),should_be)
    def assertCommand(self,should_be,index=-1):
    	self.assertEqual(self.commands[index],should_be)
    def assertCommandRegexp(self,should_be,index=-1):
    	self.assertRegexpMatches(self.commands[index],should_be)
    def test_machine_alias(self):
    	self.assertEqual(env.remote,"planck.chem.ucl.ac.uk")
    	execute(julian)
    	self.assertEqual(env.remote,"julian.chem.ucl.ac.uk")
    	execute(hector)
    	self.assertEqual(env.remote,"login.hector.ac.uk")
    def test_clean(self):
    	execute(clean)
    	self.assertCommand('make clean')
    def test_with_job(self):
        with settings(results_path="banana",local_results='pineapple'):
            with_job('foo')
            self.assertEqual(env.job_results,"banana/foo")
            self.assertEqual(env.job_results_local,"pineapple/foo")
    def test_with_template_job(self):
        with settings(results_path='banana',foo='fish',bar='swim',job_name_template="${foo}_${bar}"):     
            with_template_job()
            self.assertEqual(env.job_results,"banana/fish_swim")
    def test_hemelb(self):
        execute(hemelb,'cylinder',cores=5)
        self.assertEqual(env.name,"cylinder_abcd1234_planck_5")
        self.assertCommandRegexp('mkdir -p .*config_files/cylinder',0)
        self.assertCommandRegexp('rsync .*config_files/cylinder',1)
        self.assertCommandRegexp("put .*scripts/cylinder_abcd1234_planck_5.sh",2)
        self.assertCommandRegexp("mkdir -p .*results/cylinder_abcd1234_planck_5",3)
        self.assertCommandRegexp("cp .*scripts/cylinder_abcd1234_planck_5.sh .*results/cylinder_abcd1234_planck_5",4)
        self.assertCommandRegexp("cp .*CMakeCache.txt .*results/cylinder_abcd1234_planck_5",5)
        self.assertCommandRegexp("put .*env.yml",6)
        self.assertCommandRegexp("chmod u\+x .*scripts/cylinder_abcd1234_planck_5.sh",7)
        self.assertCommandRegexp(".*scripts/cylinder_abcd1234_planck_5.sh",8)
        self.assertCommandCount(9)
    def test_hemelbs(self):
        execute(hemelbs,'cylinder',cores='[1:6:1]')
        self.assertCommandRegexp('rsync .*config_files/cylinder',1)
        self.assertCommandRegexp("cylinder_abcd1234_planck_5.sh")
        self.assertCommandCount(9*5)
    def test_create_config(self):
        execute(create_config,'cylinder',VoxelSize=0.001)
        self.assertEqual(env.config,"cylinder_0_001_1000_3")
        self.assertCommandRegexp("mkdir -p .*/config_files/cylinder_0_001_1000_3",0)
        self.assertCommand("generate 0.001 1000 3",1)
        self.assertCommandCount(2)
    def test_create_configs(self):
        execute(create_configs,'cylinder',VoxelSize='[0.001:0.011:0.001]')
        self.assertEqual(env.config,"cylinder_0_01_1000_3")
        self.assertCommandRegexp("mkdir -p .*/config_files/cylinder_0_001_1000_3",0)
        self.assertCommand("generate 0.001 1000 3",1)
        self.assertCommandCount(2*10)
    def test_hemelb_profile(self):
        execute(hemelb_profile,'cylinder',VoxelSize='[0.001:0.011:0.001]',cores='[1:6:1]')
        self.assertEqual(env.name,"cylinder_0_01_1000_3_abcd1234_planck_5")
        self.assertCommandRegexp("mkdir -p .*/config_files/cylinder_0_001_1000_3",0)
        self.assertCommand("generate 0.001 1000 3",1)
        self.assertCommandRegexp('mkdir -p .*config_files/cylinder',2)
        self.assertCommandRegexp('rsync .*config_files/cylinder',3)
        self.assertCommandRegexp("put .*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh",4)
        self.assertCommandRegexp("mkdir -p .*results/cylinder_0_001_1000_3_abcd1234_planck_1",5)
        self.assertCommandRegexp("cp .*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh .*results/cylinder_0_001_1000_3_abcd1234_planck_1",6)
        self.assertCommandRegexp("cp .*CMakeCache.txt .*results/cylinder_0_001_1000_3_abcd1234_planck_1",7)
        self.assertCommandRegexp("put .*env.yml",8)
        self.assertCommandRegexp("chmod u\+x .*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh",9)
        self.assertCommandRegexp(".*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh",10)
        self.assertCommandCount(2*10 + 9*10*5)
    def test_hemelb_profile_no_config_generation(self):
        execute(hemelb_profile,'cylinder',VoxelSize='[0.001:0.011:0.001]',cores='[1:6:1]',create_configs="False")
        self.assertEqual(env.name,"cylinder_0_01_1000_3_abcd1234_planck_5")
        self.assertCommandRegexp('mkdir -p .*config_files/cylinder',0)
        self.assertCommandRegexp('rsync .*config_files/cylinder',1)
        self.assertCommandRegexp("put .*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh",2)
        self.assertCommandRegexp("mkdir -p .*results/cylinder_0_001_1000_3_abcd1234_planck_1",3)
        self.assertCommandRegexp("cp .*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh .*results/cylinder_0_001_1000_3_abcd1234_planck_1",4)
        self.assertCommandRegexp("cp .*CMakeCache.txt .*results/cylinder_0_001_1000_3_abcd1234_planck_1",5)
        self.assertCommandRegexp("put .*env.yml",6)
        self.assertCommandRegexp("chmod u\+x .*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh",7)
        self.assertCommandRegexp(".*scripts/cylinder_0_001_1000_3_abcd1234_planck_1.sh",8)
        self.assertCommandCount(9*10*5)
        