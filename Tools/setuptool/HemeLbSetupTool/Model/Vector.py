# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from HemeLbSetupTool.Util.Observer import Observable

class Vector(Observable):
    _Args = {'x': float("nan"),
             'y': float("nan"),
             'z': float("nan")}
    
    def __init__(self, *args):
        if len(args) == 0:
            nan = float("nan")
            self.x = nan
            self.y = nan
            self.z = nan
        elif len(args) == 1:
            self.x, self.y, self.z = args[0]
        elif len(args) == 3:
            self.x, self.y, self.z = args
        else:
            raise ValueError

        return
    
    def __str__(self):
        return '[%f, %f, %f]' % (self.x, self.y, self.z)
    def __getstate__(self):
        picdic = {}
        for attr in ('x', 'y', 'z'):
            picdic[attr] = getattr(self, attr)
            continue
        return picdic
    pass


