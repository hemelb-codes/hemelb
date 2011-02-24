class Namespace(object):
    """Simple namespace object.
    """
    def __init__(self, **kwargs):
        object.__getattribute__(self, '__dict__').update(kwargs)
        
    #     object.__setattr__(self, 'contents', {})
    #     object.__getattribute__(self, 'contents').update(kwargs)
    #     return

    def __getattribute__(self, attr):
        try:
            return object.__getattribute__(self, '__dict__')[attr]
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'" % (str(object.__getattribute__(self, '__class__')), attr))
        return

    def __setattr__(self, attr, value):
        object.__getattribute__(self, '__dict__')[attr] = value
    
    pass
