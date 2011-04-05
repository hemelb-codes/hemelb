# Site types
# cimport Flags
cdef:
    unsigned int SOLID_TYPE = 0b00
    unsigned int FLUID_TYPE = 0b01
    unsigned int INLET_TYPE = 0b10
    unsigned int OUTLET_TYPE = 0b11
    
    unsigned int BOUNDARIES = 3
    unsigned int INLET_BOUNDARY = 0
    unsigned int OUTLET_BOUNDARY = 1
    unsigned int WALL_BOUNDARY = 2

    unsigned int SITE_TYPE_BITS = 2
    unsigned int BOUNDARY_CONFIG_BITS = 14
    unsigned int BOUNDARY_DIR_BITS = 4
    unsigned int BOUNDARY_ID_BITS = 10

    unsigned int BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS
    unsigned int BOUNDARY_DIR_SHIFT = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS
    unsigned int BOUNDARY_ID_SHIFT = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS
    
    # ===============================================================================
    # Horrifying bit-fiddling masks courtesy of Marco.
    # Comments show the bit patterns.
    # ===============================================================================
    unsigned int SITE_TYPE_MASK = ((1 << SITE_TYPE_BITS) - 1)
    # 0000 0000  0000 0000  0000 0000  0000 0011
    # These give the *_TYPE given above

    unsigned int BOUNDARY_CONFIG_MASK = ((1 << BOUNDARY_CONFIG_BITS) - 1) << BOUNDARY_CONFIG_SHIFT
    # 0000 0000  0000 0000  1111 1111  1111 1100
    # These bits are set if the lattice vector they correspond to takes one to a solid site
    # The following hex digits give the index into LatticeSite.neighbours
    # ---- ----  ---- ----  DCBA 9876  5432 10--

    unsigned int BOUNDARY_DIR_MASK = ((1 << BOUNDARY_DIR_BITS) - 1) << BOUNDARY_DIR_SHIFT
    # 0000 0000  0000 1111  0000 0000  0000 0000
    # No idea what these represent. As far as I can tell, they're unused.
    
    unsigned int BOUNDARY_ID_MASK = ((1 << BOUNDARY_ID_BITS) - 1) << BOUNDARY_ID_SHIFT
    # 0011 1111  1111 0000  0000 0000  0000 0000
    # These bits together give the index of the inlet/outlet/wall in the output XML file
    
    unsigned int PRESSURE_EDGE_MASK = 1
PRESSURE_EDGE_MASK <<= 31
    # 1000 0000  0000 0000  0000 0000  0000 0000
    
