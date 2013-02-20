# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from functools import reduce
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
    
