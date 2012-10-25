#!/usr/bin/env python
# encoding: utf-8
"""
test_utils.py

Created by James Hetherington on 2012-10-18.
Copyright (c) 2012 UCL. All rights reserved.
"""

import sys
import os
import unittest
import numpy as np
from ..utils import *

class utilsTests(unittest.TestCase):
	def setUp(self):
		pass
	def testMatchCorresponding(self):
	    x=np.array([1,2,7,3,28])
	    y=np.array([2,28,7,3,1])
	    result=MatchCorresponding(x,y)
	    self.assertSequenceEqual(list(result),[1,4,2,3,0])

if __name__ == '__main__':
	unittest.main()
