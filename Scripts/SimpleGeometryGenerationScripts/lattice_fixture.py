# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import zlib
import xdrlib
from functools import reduce
from Site import Site
from link import Link

class LatticeFixture(object):

    # Versioning constants
    hemelb_magic_number = 0x686c6221
    geometry_magic_number = 0x676d7904
    hemelb_geometry_file_format_version_number = 4

    def __init__(self):
        # Assuming zero lattice site is at the origin of coordinates (0,0,0)
        self.origin = (0.0,) * 3
        self.pack = xdrlib.Packer()

    def write_preamble(self):
        preamble_encoder = xdrlib.Packer()
        # Magic numbers and versioning
        preamble_encoder.pack_uint(self.hemelb_magic_number) 
        preamble_encoder.pack_uint(self.geometry_magic_number)
        preamble_encoder.pack_uint(self.hemelb_geometry_file_format_version_number)

        # Domain size in blocks
        for blocks_across_dimension in self.size_in_blocks:
            preamble_encoder.pack_uint(blocks_across_dimension)

        # Number of lattice sites along the side of a block
        preamble_encoder.pack_uint(self.sites_along_block)

        # Pad the length of the preamble to 64 bytes
        preamble_encoder.pack_uint(0)
        self.outfile.write(preamble_encoder.get_buffer())

    def write_dummy_header(self):
        header_encoder = xdrlib.Packer()
        for block in self.blocks:
            header_encoder.pack_uint(0)
            header_encoder.pack_uint(0)
            header_encoder.pack_uint(0)
            continue

        self._header_start = self.outfile.tell()
        self.outfile.write(header_encoder.get_buffer())
        self._header_end = self.outfile.tell()
        return

    def rewrite_header(self):
        header_encoder = xdrlib.Packer()
        for block in self.blocks:
            header_encoder.pack_uint(block.get_number_of_fluid_sites())
            header_encoder.pack_uint(block.compressed_bytes)
            header_encoder.pack_uint(block.uncompressed_bytes)
            continue

        self.outfile.seek(self._header_start)
        self.outfile.write(header_encoder.get_buffer())
        assert self.outfile.tell() == self._header_end
        return

    def write_blocks(self):
        for block in self.blocks:
            if block.get_number_of_fluid_sites() == 0:
                block.compressed_bytes = 0
                block.uncompressed_bytes = 0
                continue

            block_encoder = xdrlib.Packer()
            for site in block.sites:
                self.pack_site(site, block_encoder)
                continue

            compressed = zlib.compress(block_encoder.get_buffer())
            block.compressed_bytes = len(compressed)
            block.uncompressed_bytes = len(block_encoder.get_buffer())

            self.outfile.write(compressed)
            continue
        return

    def pack_site(self, site, block_encoder):
        block_encoder.pack_uint(site.site_type)
        # If the site is solid, nothing else is packed. Otherwise, pack information about the links
        if site.site_type == Site.fluid_site:
            for link in site.links:
                self.pack_link(link, block_encoder)
            if site.normal is None:
                block_encoder.pack_uint(0)
            else:
                block_encoder.pack_uint(1)
                block_encoder.pack_float(site.normal[0])
                block_encoder.pack_float(site.normal[1])
                block_encoder.pack_float(site.normal[2])

    def pack_link(self, link, block_encoder):
        # Pack the link type and depending on it, extra information
        block_encoder.pack_uint(link.link_type)
        if link.link_type in [Link.inlet, Link.outlet]:
            # An unsigned integer giving the an index to the inlet array specified in the XML Config File
            block_encoder.pack_uint(link.iolet_index)
        if link.link_type != Link.no_boundary:
            # A single precision floating point number giving the distance to the boundary, as a fraction of the lattice vector
            block_encoder.pack_float(link.wall_distance)

    def write(self, filename):
        self.outfile = file(filename, 'wb')
        self.write_preamble()
        self.write_dummy_header()
        self.write_blocks()
        self.rewrite_header()
        self.outfile.close()
        return

    lattice_directions = (
        (-1, -1, -1),
        (-1, -1, 0),
        (-1, -1, +1),
        (-1, 0, -1),
        (-1, 0, 0),
        (-1, 0, +1),
        (-1, +1, -1),
        (-1, +1, 0),
        (-1, +1, +1),
        (0, -1, -1),
        (0, -1, 0),
        (0, -1, +1),
        (0, 0, -1),
        #( 0, 0, 0),
        (0, 0, +1),
        (0, +1, -1),
        (0, +1, 0),
        (0, +1, +1),
        (+1, -1, -1),
        (+1, -1, 0),
        (+1, -1, +1),
        (+1, 0, -1),
        (+1, 0, 0),
        (+1, 0, +1),
        (+1, +1, -1),
        (+1, +1, 0),
        (+1, +1, +1)
    )
