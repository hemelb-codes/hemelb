#!/usr/bin/env python
# encoding: utf-8
"""
test_machine_environment.py

Created by James Hetherington on 2012-01-19.
Copyright (c) 2012 UCL. All rights reserved.
"""
import unittest
import sys

from ..fab import *

class TestFabric(unittest.TestCase):
	def setUp(self):
		#Update the user config with testing example
		user_config=yaml.load(open(os.path.join(env.localroot,'deploy','test','machines_user.yml')))
		execute(planck) #Default machine target is assumed as planck.
		#Monkeypatch the fabric run command to push the command onto a mock command list.
		sys.modules['deploy.fab'].run=lambda command: self.commands.append(command)
		self.commands=[]
	def assertLastCommand(self,should_be):
		self.assertEqual(self.commands[-1],should_be)
	def test_machine_alias(self):
		self.assertEqual(env.remote,"planck.chem.ucl.ac.uk")
		execute(julian)
		self.assertEqual(env.remote,"julian.chem.ucl.ac.uk")
		execute(hector)
		self.assertEqual(env.remote,"login.hector.ac.uk")
	def test_clean(self):
		execute(clean)
		self.assertLastCommand('make clean')
		
