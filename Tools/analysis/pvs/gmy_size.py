import numpy as np
from hemeTools.parsers.geometry import ConfigLoader

class GotSizeException(Exception):
    def __init__(self, size):
        self.size = size
        return
    pass
    
class SizeGetter(ConfigLoader):
    def OnEndHeader(self):
        raise GotSizeException(np.sum(self.Domain.BlockFluidSiteCounts))
    pass

def GetGmySize(gmy):
    getter = SizeGetter(gmy)
    try:
        getter.Load()
    except GotSizeException as e:
        return e.size
    else:
        raise RuntimeError("SizeGetter should have thrown a GotSizeException!")
    
if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        print GetGmySize(arg)
