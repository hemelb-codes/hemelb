
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
"""Memoize computed properties.

"""

class memo_property(object):
    """Decorator for memoizing computed, static properties.
    """
    def __init__(self, fget):
        self.memo_attr_name = '_memo_' + fget.func_name
        self.getter = fget
        return
    
    def __get__(self, instance, owner_cls):
        try:
            ans = getattr(instance, self.memo_attr_name)
        except AttributeError:
            ans = self.getter(instance)
            setattr(instance, self.memo_attr_name, ans)
            pass
        return ans
    pass
