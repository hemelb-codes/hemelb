# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import yaml

class LoadableMixin:
    '''Define a class method from_dict and this will add loading from
    streams and files.
    '''
    
    @classmethod
    def from_file(cls, fn):
        with open(fn) as f:
            return cls.from_stream(f)
        
    @classmethod
    def from_stream(cls, stream):
        state = yaml.load(stream)
        return cls.from_dict(state)

    pass

class DumpableMixin:
    '''Define a method to_dict and this will add writing to streams and
    files.
    '''
    def to_file(self, fn):
        with open(fn, 'w') as f:
            self.to_stream(f)

    def to_stream(self, stream):
        state = self.to_dict()
        yaml.dump(state, stream)
        
    pass

__all__ = ['LoadableMixin', 'DumpableMixin']
