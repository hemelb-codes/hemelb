#!/usr/bin/env python
# encoding: utf-8
"""
test_paged_socket.py

Created by James Hetherington on 2012-04-16.
Copyright (c) 2012 University College London. All rights reserved.
"""

import unittest
import mock

from hemelb_steering.paged_socket import PagedSocket

class TestPagedSocket(unittest.TestCase):
    def setUp(self):
        with mock.patch('socket.socket') as mock_socket_maker:
            self.mockSocket = mock_socket_maker.return_value
            self.ps = PagedSocket(address='fibble', port=8080, receive_length=9)
            mock_socket_maker.assert_called_once_with()
            self.mockSocket.connect.assert_called_once_with(('fibble', 8080))
            
    def test_connect(self):
        pass #The mocks defined in setUp are sufficient here
        
    def test_receive(self):
        self.mockSocket.recv.return_value = 'dummyPage'
        self.assertEqual(bytearray('dummyPage'), self.ps.receive())
        self.mockSocket.recv.assert_called_once_with(9)
        
    def test_send(self):
        payload = str(bytearray([1, 2, 3, 4, 5]))
        self.mockSocket.send.return_value = len(payload)
        self.ps.send(payload)
        self.mockSocket.send.assert_called_once_with(payload)
        
    def test_staged_send(self):
        payload=str(bytearray([1, 2, 3, 4, 5]))
        self.mockSocket.send.side_effect = [3, 2]
        self.ps.send(payload)
        self.mockSocket.send.assert_has_calls([mock.call(str(bytearray([1, 2, 3, 4, 5]))), mock.call(str(bytearray([4, 5])))])
        
    def test_variable_length_receive(self):
    
        def _dummy_additional_receive(header):
            return bytearray(header)[0]
            
        self.ps.additional_receive_length_function = _dummy_additional_receive
        self.mockSocket.recv.side_effect = [str(bytearray([5, 2, 3, 4, 5, 6, 7, 8, 9])), str(bytearray([10, 11, 12, 13, 14]))]
        self.assertEqual(bytearray([5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]), self.ps.receive())
        self.mockSocket.recv.assert_has_calls([mock.call(9), mock.call(5)])
        
        
        
