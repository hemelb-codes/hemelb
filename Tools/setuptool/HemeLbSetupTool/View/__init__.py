from MainWindow import MainWindow
from Delegator import Delegator

class Namespace(object):

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

class WxView(Delegator):

    def __init__(self):
        self.exports = Namespace()
        Delegator.__init__(self, self.exports)
        return

    def Create(self, controller):
        self.controller = controller
        self.main = MainWindow(self)
        
        return self.main
    
