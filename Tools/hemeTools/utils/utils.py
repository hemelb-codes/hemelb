#!/usr/bin/env python
# encoding: utf-8
"""
utils.py

Created by James Hetherington on 2012-10-18.
Copyright (c) 2012 UCL. All rights reserved.

Utility functions for HemeTools
"""

import numpy as np

def MatchCorresponding(first, second):
    """
    Given two 1D arrays of the same length, containing different
    permutations of the same numbers, return an array results such
    that first[result[i]] == second[i]
    """
    assert first.shape == second.shape
    assert first.ndim == 1
    n = len(first)
    
    fAS = np.argsort(first)
    sAS = np.argsort(second)
    # Know that first[fAS[i]] == second[sAS[i]]
    result = np.empty_like(second)
    for i in xrange(n):
        result[sAS[i]] = fAS[i]
    
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
