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
