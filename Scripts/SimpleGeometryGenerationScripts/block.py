# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from Site import Site

class Block(object):
    def __init__(self, sites):
        self.sites = sites

    def get_number_of_fluid_sites(self):
        n = 0
        for site in self.sites:
            if site.site_type == Site.fluid_site:
                n += 1
        return n
    
