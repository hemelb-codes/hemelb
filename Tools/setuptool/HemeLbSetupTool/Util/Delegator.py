
class Delegator(object):
    """Mixin that makes the object delegate unknown attributes to the
    specified delegate.
    """
    def __init__(self, delegate):
        self._delegate = delegate
        return
    
    def __getattr__(self, attr):
        """Delegate unknown attributes to the parent.
        """
        try:
            return getattr(object.__getattribute__(self, '_delegate'),
                           attr)
        except AttributeError:
            raise AttributeError("'%s' object (nor its parents) has no attribute '%s'" % (str(object.__getattribute__(self, '__class__')), attr))
        return
    pass


class ParentDelegator(Delegator):
    """Mixin that makes the object delegate to its parents.
    """
    def __init__(self):
        self._delegate = self.GetParent()
        return
