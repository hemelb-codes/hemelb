#!/usr/bin/env python
# encoding: utf-8
"""
utils.py

Created by James Hetherington on 2012-10-18.
Copyright (c) 2012 UCL. All rights reserved.

Utility functions for HemeTools
"""

import sys
import os
import numpy as np

def MatchCorresponding(first,second):
    count=second.shape[0]
    result=np.array([np.argwhere(first==val)[0] for val in second])
    return result

def MatchCorrespondingOld(first,second):
    count=second.shape[0]
    result = np.empty_like(first)
    
    for index in xrange(count):
        for search in xrange(count):
            if first[search] == second[index]:
                result[index] = search
                break
                
    return result