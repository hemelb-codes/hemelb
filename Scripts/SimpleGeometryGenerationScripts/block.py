from functools import reduce
from Site import Site
class Block(object):
    def __init__(self,sites):
        self.sites=sites
        self.bytes=sum(map(Site.bytes,self.sites))
