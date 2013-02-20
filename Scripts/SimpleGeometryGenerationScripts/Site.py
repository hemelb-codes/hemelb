# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

class Site(object):
    # Type of site identifier
    solid_site, fluid_site = range(2)

    # Bytes stored per link type:
    #   none:   link type (uint)
    #   wall:   link type and distance to boundary (uint+float),
    #   inlet:  link type and inlet index and distance to boundary (uint+uint+float)
    #   outlet: link type and inlet index and distance to boundary (uint+uint+float)
    link_extra_bytes_required = [4, 4 + 4, 4 + 4 + 4, 4 + 4 + 4]
    
    def __init__(self, site_type, links):
        self.site_type = site_type
        self.links = links

    def bytes(self):
        # unsigned site identifier (solid vs fluid) plus the data required to represent each link
        size = 4 + sum(self.link_extra_bytes_required[link.link_type] for link in self.links)
        return size
