# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

class Site(object):
    # Type of site identifier
    solid_site, fluid_site = range(2)

    # Bytes stored per link type:
    #   none:   link type (uint)
    #   wall:   link type and distance to boundary (uint+float),
    #   inlet:  link type and inlet index and distance to boundary (uint+uint+float)
    #   outlet: link type and inlet index and distance to boundary (uint+uint+float)
    link_extra_bytes_required = [4, 4 + 4, 4 + 4 + 4, 4 + 4 + 4]
    
    def __init__(self, site_type, links, normal):
        self.site_type = site_type
        self.links = links
        self.normal = normal

    def bytes(self):
        # unsigned site identifier (solid vs fluid) plus the data required to represent each link
        # plus an unsigned telling whether there's a normal available
        size = 4
        if self.site_type == self.fluid_site:
            size += sum(self.link_extra_bytes_required[link.link_type] for link in self.links) + 4
        
            # if there's a wall normal available, we need 3 floats to store them
            if self.normal is not None:
                size += 3 * 4
            
        return size
