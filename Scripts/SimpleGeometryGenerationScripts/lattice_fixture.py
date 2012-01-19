import xdrlib
from functools import reduce
from Site import Site
from link import Link

class LatticeFixture(object):

    # Versioning constants
    hemelb_magic_number = 0x686c6221
    geometry_magic_number = 0x676d7904
    hemelb_geometry_file_format_version_number = 2

    def __init__(self):
        # Assuming zero lattice site is at the origin of coordinates (0,0,0)
        self.origin=(0.0,)*3
        self.pack=xdrlib.Packer()

    def pack_all(self):
        self.pack_preamble()
        self.pack_header()
        self.pack_blocks()
        
    def pack_preamble(self):
        # Magic numbers and versioning
        self.pack.pack_uint(self.hemelb_magic_number) 
        self.pack.pack_uint(self.geometry_magic_number)
        self.pack.pack_uint(self.hemelb_geometry_file_format_version_number)

        # Domain size in blocks
        for blocks_across_dimension in self.size_in_blocks:
            self.pack.pack_uint(blocks_across_dimension)

        # Number of lattice sites along the side of a block
        self.pack.pack_uint(self.sites_along_block)

        # Grid space step
        self.pack.pack_double(self.space_step)

        # Origin of coordinates
        for origin_coordinate in self.origin:
            self.pack.pack_double(origin_coordinate)

        # Pad the length of the preamble to 64 bytes
        self.pack.pack_uint(0)

    def pack_header(self):
        for block in self.blocks:
            self.pack_block_header(block)

    def pack_block_header(self,block):
        self.pack.pack_uint(len(block.sites))
        self.pack.pack_uint(block.bytes)

    def pack_blocks(self):
        # Requirement: TODO Data is compressed with $COMPRESSION_ALGORITHM on a per-block basis.
        for block in self.blocks:
            start_position = len(self.pack.get_buffer())
            # Requirement: Only if a block contains at least one fluid site is data given for the block
            for site in block.sites:
                self.pack_site(site)
            end_position = len(self.pack.get_buffer())
            assert end_position-start_position == block.bytes, "Inconsistent block size evaluation"

    def pack_site(self,site):
        self.pack.pack_uint(site.site_type)
        # If the site is solid, nothing else is packed. Otherwise, pack information about the links
        if site.site_type == Site.fluid_site:
            for link in site.links:
                self.pack_link(link)

    def pack_link(self,link):
        # Pack the link type and depending on it, extra information
        self.pack.pack_uint(link.link_type)
        if link.link_type in [Link.inlet,Link.outlet]:
            # An unsigned integer giving the an index to the inlet array specified in the XML Config File
            self.pack.pack_uint(link.iolet_index)
        if link.link_type != Link.no_boundary:
            # A single precision floating point number giving the distance to the boundary, as a fraction of the lattice vector
            self.pack.pack_float(link.wall_distance)

    def write(self,filename):
        self.pack_all()
        file_handle = file(filename,'wb')
        file_handle.write(self.pack.get_buffer())

    lattice_directions=(
        (-1,-1,-1),
        (-1,-1, 0),
        (-1,-1,+1),
        (-1, 0,-1),
        (-1, 0, 0),
        (-1, 0,+1),
        (-1,+1,-1),
        (-1,+1, 0),
        (-1,+1,+1),
        ( 0,-1,-1),
        ( 0,-1, 0),
        ( 0,-1,+1),
        ( 0, 0,-1),
        #( 0, 0, 0),
        ( 0, 0,+1),
        ( 0,+1,-1),
        ( 0,+1, 0),
        ( 0,+1,+1),
        (+1,-1,-1),
        (+1,-1, 0),
        (+1,-1,+1),
        (+1, 0,-1),
        (+1, 0, 0),
        (+1, 0,+1),
        (+1,+1,-1),
        (+1,+1, 0),
        (+1,+1,+1)
    )
