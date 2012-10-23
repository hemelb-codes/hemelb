#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 


from lattice_fixture import LatticeFixture
from link import Link
from Site import Site
from block import Block
import numpy as np

class FourCube(LatticeFixture):
    def __init__(self):
        # Domain is (0,3)x(0,3)X(0,3) length unit
        x_min,y_min,z_min=0,0,0
        x_max,y_max,z_max=3,3,3
        
        # Single block with 4x4x4 sites, 1 length unit apart
        self.size_in_blocks=(1,1,1)
        self.sites_along_block=4
        self.space_step=1

        sites=[]
        # Requirement: Sites are striped with the z coordinate changing most frequently, i.e. the z coordinate represents the least significant part of the site index within the block.
        for index in np.ndindex(3*(self.sites_along_block,)):
            x_index,y_index,z_index=index
            links=[]
            for direction in self.lattice_directions:
                link_type = Link.no_boundary

                # Assign an iolet id (the first iolet in the config xml file)
                iolet_index=0

                # In this example, links that cross both wall and inlet/outlet (like the 8 corners of the cube) will be considered non-wall because of the ordering of the if statements below        
                if x_index==x_min and direction[0]==-1:
                    link_type = Link.wall
                if x_index==x_max and direction[0]==1:
                    link_type = Link.wall
                if y_index==y_min and direction[1]==-1:
                    link_type = Link.wall
                if y_index==y_max and direction[1]==1:
                    link_type = Link.wall
                if z_index==z_min and direction[2]==-1:
                    link_type = Link.inlet
                if z_index==z_max and direction[2]==1:
                    link_type = Link.outlet

                # Assume walls are half a lattice away
                wall_distance=0.5

                # iolet_index and wall_distance will be ignored when meaningless
                links.append(Link(link_type,wall_distance,iolet_index))

            sites.append(Site(Site.fluid_site,links))

        # Requirement: Blocks are striped with the z coordinate changing most frequently (i.e. the z coordinate represents the least significant part of the block index) then y and x slowest.
        # This geometry is made of a single block with the sites defined above
        self.blocks=[Block(sites)]

        # Let the parent class do the rest
        super(FourCube,self).__init__()

if __name__ == '__main__':
    # TODO: get the filename from command line
    filename = "four_cube.gmy"
    cube = FourCube()
    cube.write(filename)