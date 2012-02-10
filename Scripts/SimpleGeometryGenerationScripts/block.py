from functools import reduce
from Site import Site

class Block(object):
    def __init__(self,sites):
        self.sites=sites

    def get_number_of_fluid_sites(self):
        n = 0
        for site in self.sites:
            if site.site_type == Site.fluid_site:
                n += 1
        return n
    
