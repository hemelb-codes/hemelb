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
import textwrap

from ..fab import *

class TestFabric(unittest.TestCase):
    def setUp(self):
    	#Update the user config with testing example
    	env.test_home=os.path.join(env.localroot,'deploy','test')
    	user_config=yaml.load(open(os.path.join(env.localroot,'deploy','test','machines_user.yml')))
    	env.update(user_config['default'])
    	execute(planck) #Default machine target is assumed as planck.
    	#Monkeypatch the fabric commands to do nothing, but record what they would have done
    	sys.modules['deploy.fab'].run=lambda command: self.commands.append(command)
    	def mock_local(command,original=sys.modules['deploy.fab'].local):
    	  self.commands.append(command)
    	  original(command)
    	sys.modules['deploy.fab'].local=mock_local  
    	sys.modules['deploy.fab'].put=lambda source,target: self.commands.append("put "+source+" "+target)
    	sys.modules['deploy.fab'].rsync_project=lambda **args: self.commands.append("rsync "+args['local_dir']+" "+args['remote_dir'])
    	def mock_profile(profile,original=sys.modules['deploy.fab'].generate):
    	   self.commands.append("generate %g %g %g"%(profile.VoxelSize, profile.Steps , profile.Cycles) )
    	   original(profile)
    	sys.modules['deploy.fab'].generate=mock_profile
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
        self.assertEqual(env.name,"cylinder_abcd1234_planck_5_10_10")
        self.assertCommandRegexp('mkdir -p .*config_files/cylinder',0)
        self.assertCommandRegexp('rsync .*config_files/cylinder',1)
        self.assertCommandRegexp("put .*scripts/cylinder_abcd1234_planck_5_10_10.sh",2)
        self.assertCommandRegexp("mkdir -p .*results/cylinder_abcd1234_planck_5_10_10",3)
        self.assertCommandRegexp("cp .*scripts/cylinder_abcd1234_planck_5_10_10.sh .*results/cylinder_abcd1234_planck_5_10_10",4)
        self.assertCommandRegexp("cp .*CMakeCache.txt .*results/cylinder_abcd1234_planck_5_10_10",5)
        self.assertCommandRegexp("put .*env.yml",6)
        self.assertCommandRegexp("chmod u\+x .*scripts/cylinder_abcd1234_planck_5_10_10.sh",7)
        self.assertCommandRegexp(".*scripts/cylinder_abcd1234_planck_5_10_10.sh",8)
        self.assertCommandCount(9)
    def test_hemelbs(self):
        execute(hemelbs,'cylinder',cores='[1:6:1]')
        self.assertCommandRegexp('rsync .*config_files/cylinder',1)
        self.assertCommandRegexp("cylinder_abcd1234_planck_5_10_10.sh")
        self.assertCommandCount(9*5)
    def test_create_config(self):
        execute(create_config,'cylinder',VoxelSize=0.1)
        self.assertEqual(env.config,"cylinder_0_1_1000_3")
        self.assertCommandRegexp("mkdir -p .*/configs/cylinder_0_1_1000_3",0)
        self.assertCommand("generate 0.1 1000 3",1)
        self.assertCommandCount(2)
    def test_create_configs(self):
        execute(create_configs,'cylinder',VoxelSize='[0.1:0.21:0.01]')
        self.assertEqual(env.config,"cylinder_0_2_1000_3")
        self.assertCommandRegexp("mkdir -p .*/configs/cylinder_0_1_1000_3",0)
        self.assertCommand("generate 0.1 1000 3",1)
        self.assertCommandCount(2*11)
    def test_hemelb_profile(self):
        execute(hemelb_profile,'cylinder',VoxelSize='[0.1:0.21:0.01]',cores='[1:6:1]')
        self.assertEqual(env.name,"cylinder_0_2_1000_3_abcd1234_planck_5_10_10")
        self.assertCommandRegexp("mkdir -p .*/configs/cylinder_0_1_1000_3",0)
        self.assertCommand("generate 0.1 1000 3",1)
        self.assertCommandRegexp('mkdir -p .*config_files/cylinder',2)
        self.assertCommandRegexp('rsync .*config_files/cylinder',3)
        self.assertCommandRegexp("put .*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh",4)
        self.assertCommandRegexp("mkdir -p .*results/cylinder_0_1_1000_3_abcd1234_planck_1_10_10",5)
        self.assertCommandRegexp("cp .*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh .*results/cylinder_0_1_1000_3_abcd1234_planck_1_10_10",6)
        self.assertCommandRegexp("cp .*CMakeCache.txt .*results/cylinder_0_1_1000_3_abcd1234_planck_1_10_10",7)
        self.assertCommandRegexp("put .*env.yml",8)
        self.assertCommandRegexp("chmod u\+x .*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh",9)
        self.assertCommandRegexp(".*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh",10)
        self.assertCommandCount(2*11 + 9*11*5)
    def test_hemelb_profile_no_config_generation(self):
        execute(hemelb_profile,'cylinder',VoxelSize='[0.1:0.21:0.01]',cores='[1:6:1]',create_configs="False")
        self.assertEqual(env.name,"cylinder_0_2_1000_3_abcd1234_planck_5_10_10")
        self.assertCommandRegexp('mkdir -p .*config_files/cylinder',0)
        self.assertCommandRegexp('rsync .*config_files/cylinder',1)
        self.assertCommandRegexp("put .*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh",2)
        self.assertCommandRegexp("mkdir -p .*results/cylinder_0_1_1000_3_abcd1234_planck_1_10_10",3)
        self.assertCommandRegexp("cp .*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh .*results/cylinder_0_1_1000_3_abcd1234_planck_1_10_10",4)
        self.assertCommandRegexp("cp .*CMakeCache.txt .*results/cylinder_0_1_1000_3_abcd1234_planck_1_10_10",5)
        self.assertCommandRegexp("put .*env.yml",6)
        self.assertCommandRegexp("chmod u\+x .*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh",7)
        self.assertCommandRegexp(".*scripts/cylinder_0_1_1000_3_abcd1234_planck_1_10_10.sh",8)
        self.assertCommandCount(9*11*5)
    def test_configure_default(self):
        execute(configure)
        target={
            'CMAKE_BUILD_TYPE': "Release",
            'CMAKE_CXX_FLAGS_RELEASE': "-O4",
            'CMAKE_INSTALL_PREFIX': env.install_path,
            'CPPUNIT_PATCH_LDL' : True,
            "HEMELB_DEPENDENCIES_INSTALL_PATH": env.install_path,
            "HEMELB_SUBPROJECT_MAKE_JOBS": 1
        }
        self.assertEqual(env.total_cmake_options,target)
        #Can't just assert on a string here, as the order of the dict is not defined
        for key,value in target.iteritems():
            self.assertRegexpMatches(env.cmake_flags,"-D%s=%s"%(key,value))
    def test_configure_debug(self):
        execute(configure,'debug')
        self.assertEqual(env.total_cmake_options,
        {
            'CMAKE_BUILD_TYPE': "Debug",
            'HEMELB_OPTIMISATION': "",
            'HEMELB_LOG_LEVEL': "debug",
            'CPPUNIT_PATCH_LDL' : True,
            'CMAKE_INSTALL_PREFIX': env.install_path,
            "HEMELB_DEPENDENCIES_INSTALL_PATH": env.install_path,
            "HEMELB_SUBPROJECT_MAKE_JOBS": 1
        })
        
    def test_script_template(self):
      script=script_templates('dummy_ge_header','dummy_jobscript',commands=['extra'])
      content=open(script).read()
      self.assertEqual(content,"user: test_user\n\nrun bananas\n\nextra")
