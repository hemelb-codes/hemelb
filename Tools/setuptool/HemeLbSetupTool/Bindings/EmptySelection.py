# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from ..Util.Observer import Observable

class EmptySelection(Observable):
    """Represent an empty selection.
    """
    def __getattr__(self, attr):
        return self
    
    pass

EmptySelection = EmptySelection()


def isNone(value):
    return (value is None) or (value is EmptySelection)
