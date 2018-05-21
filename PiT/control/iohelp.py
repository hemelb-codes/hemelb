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
