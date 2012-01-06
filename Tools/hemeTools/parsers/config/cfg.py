# ======================================================================
# Bit fiddling helpers for the "cfg" variable in the files
# ======================================================================
def GetType(cfg):
    return cfg & SITE_TYPE_MASK
def GetBoundaryConfig(cfg):
    return (cfg & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT
def GetBoundaryDir(cfg):
    return (cfg & BOUNDARY_DIR_MASK) >> BOUNDARY_DIR_SHIFT
def GetBoundaryId(cfg):
    return (cfg & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT
def GetPressureEdge(cfg):
    return cfg & PRESSURE_EDGE_MASK
SOLID_TYPE = 0b00
FLUID_TYPE = 0b01
INLET_TYPE = 0b10
OUTLET_TYPE = 0b11

BOUNDARIES = 3
INLET_BOUNDARY = 0
OUTLET_BOUNDARY = 1
WALL_BOUNDARY = 2

SITE_TYPE_BITS = 2
BOUNDARY_CONFIG_BITS = 26
BOUNDARY_DIR_BITS = 4
BOUNDARY_ID_BITS = 10

BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS
BOUNDARY_DIR_SHIFT = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS
BOUNDARY_ID_SHIFT = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS

# Comments show the bit patterns of the preceeding variable

SITE_TYPE_MASK = ((1 << SITE_TYPE_BITS) - 1)
# 0000 0000  0000 0000  0000 0000  0000 0011
# These give the *_TYPE given above

BOUNDARY_CONFIG_MASK = ((1 << BOUNDARY_CONFIG_BITS) - 1) << BOUNDARY_CONFIG_SHIFT
# 0000 0000  0000 0000  1111 1111  1111 1100
# These bits are set if the lattice vector they correspond to takes one to a solid site
# The following hex digits give the index into LatticeSite.neighbours
# ---- ----  ---- ----  DCBA 9876  5432 10--

BOUNDARY_DIR_MASK = ((1 << BOUNDARY_DIR_BITS) - 1) << BOUNDARY_DIR_SHIFT
# 0000 0000  0000 1111  0000 0000  0000 0000
# No idea what these represent. As far as I can tell, they're unused.

BOUNDARY_ID_MASK = ((1 << BOUNDARY_ID_BITS) - 1) << BOUNDARY_ID_SHIFT
# 0011 1111  1111 0000  0000 0000  0000 0000
# These bits together give the index of the inlet/outlet/wall in the output XML file

PRESSURE_EDGE_MASK = 1 << (BOUNDARY_ID_SHIFT + BOUNDARY_ID_BITS + 1)
# 1000 0000  0000 0000  0000 0000  0000 0000
