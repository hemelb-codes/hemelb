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
	    x=np.array([1,2,3,3,5])
	    y=np.array([2,3,5,5,1,])
	    result=MatchCorrespondingB(x,y)
	    self.assertSequenceEqual(list(result),[1,2,4,4,0])

if __name__ == '__main__':
	unittest.main()