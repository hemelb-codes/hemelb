#!/usr/bin/env python
# encoding: utf-8
"""
test_remote_hemelb.py

Created by James Hetherington on 2012-04-16.
Copyright (c) 2012 University College London. All rights reserved.
"""

import unittest
import mock

from hemelb_steering.remote_hemelb import RemoteHemeLB 


class TestRemoteHemeLB(unittest.TestCase):
	def setUp(self):
	    with mock.patch('hemelb_steering.remote_hemelb.PagedSocket') as mockPagedSocket:
	        self.mockSocket=mockPagedSocket.return_value
	        self.rhlb=RemoteHemeLB(address='fibble',port=8080,steering_id=1111)     
	        mockPagedSocket.assert_called_once_with(address='fibble',port=8080)
	def test_receive_image(self):
	    self.mockSocket.receive.return_value='dummyPage'
	    self.rhlb.step()
	    self.assertEqual('dummyPage',self.rhlb.image)
	def test_two_step_one_image(self):
	    self.mockSocket.receive.return_value='dummyPage'
	    self.rhlb.step()
	    self.assertEqual('dummyPage',self.rhlb.image)
	    self.mockSocket.receive.assert_called_once_with()
	    self.mockSocket.receive.return_value=None
	    self.rhlb.step()
	    self.assertEqual('dummyPage',self.rhlb.image)
	    self.mockSocket.receive.assert_has_calls([mock.call(),mock.call()])
	def test_no_change_no_send(self):
	    self.rhlb.step()
	    self.assertEqual(self.mockSocket.send.mock_calls,[])
	def test_change_send(self):
	    self.rhlb.latitude=50
	    self.rhlb.step() 
	    self.mockSocket.send.assert_called_once_with(50)