# Ancient outdated geometry format

This format is not used, but we're keeping this in case we need to recover one of theses.

Please see the [new format description](geometry.md) instead.

## Geometry File Format

This page describes the XDR file format used to describe HemeLB lattice geometries.

### Description

The file is binary data using the [XDR standard](http://tools.ietf.org/html/rfc4506). This is read/written using standard libraries, in both C, using [<rpc/xdr.h>](http://linux.die.net/man/3/xdr), and in Python, using [xdrlib](http://docs.python.org/library/xdrlib.html).

This file describes the problem domain as a series of  blocks, each block being an identically sized cubic subsection of the problem domain.  

#### Preamble


The file begins with a preamble. This consists of:

* An unsigned integer representing the stress type HemeLB should use when calculating stress.
  * 0 Von Mises Stress
  * 1 Shear Stress
  * 2 Ignore Stress
* Three unsigned integers representing the x, y, and z size of the problem domain, in blocks.
* An unsigned integer representing the number of lattice sites along the side of a block. (So the total sites in a block is the cube of this.)
* A double precision floating point number representing the side length of a lattice site (or voxel), in the problem domain coordinate system.
* Three double precision floating point numbers representing the x,y, and z coordinates, in the problem domain coordinate system, of the zero lattice site.

#### Block Headers

* The file then contains an array of pairs of unsigned integers, describing each block. 
* Blocks are striped with the z coordinate changing most frequently, i.e. the z coordinate represents the least significant part of the block index.

Each pair is:

* An unsigned integer giving the number of fluid lattice sites in the block.
  * This does not have to be the same as the cube of the side length of a block in lattice sites, because some of the sites in the block might be "solid" sites which are not simulated in the Lattice Boltzmann.
* An unsigned integer giving the total length, in bytes, of the data representing a given block
  * This is useful to allow one to jump to the data for a given block (by computing the partial sum of these offsets up to the block index).

#### Block data

* Only if a block contains at least one fluid site is data given for the block. 
* If a block has no fluid sites, both the number of sites in the block, and the length of the data for the block, as given in the block headers array, will be zero. 
* A parser must use this knowledge so that the index of the block for a given sequence of block data is established correctly.

* The data for a block consists of data for each site in the block. 
* Sites are striped with the z coordinate changing most frequently, i.e. the z coordinate represents the least significant part of the site index within the block. 
* For each site, the data given, and the length of the data given, depends on the type of the site. 
* The data is given site by site, one after the other, without padding. Because of this, it is not possible to find the location in the block data of the data for a particular site without reading all the preceding sites.

#### Sites

For each site type, the first data given is an unsigned integer giving the type of the block. 
This can be used to deduce the length of the site data, and the start of the location of the next site. 
This is a bitwise structure describing several things about a site.

* The two least significant bits describe 
  * 0 Solid site
  * 1 Fluid site
  * 2 Inlet site
  * 3 Outlet site
* The next fourteen bits describe whether there is a wall or boundary in each of the fourteen non-stationary directions of the D3Q15 Lattice Boltzmann lattice.
  * The bit is set to one if there is a boundary
  * The directions are ordered as in source:Code/D3Q15.cc (the setuptool links against this code)
* The next four bits are unused. They _may_ have previously encoded the direction closest to the boundary normal.
* The next ten bits give an index into either the inlet array or outlet array specified in the [XML Config File](XMLConfigFile/), depending on whether the site is an inlet or outlet.
  * This puts an upper limit of 1024 on the number of inlets and outlets which may exist.
* The final bit specifies whether the site is a pressure edge, i.e. any of the directions intersect a wall

If the site is a simple fluid or solid site, i.e. it is not near a wall or inlet or outlet (iolet) then this completes the site data.

If the site is an inlet or outlet (bit two of the configuration byte is set) then there follows data describing the boundary
* Three double precision floating point values describing the normal to the boundary
  * The sense of the normal is given out of the fluid  into the solid
* A double precision floating point value giving the distance to the boundary, in lattice units.

If the site is adjacent to a wall then there follows data describing the wall
* Three double precision floating point values describing the normal to wall (at the closest point of intersection of the lattice vectors with the wall)
  * The sense of the normal is given out of the fluid into the solid
* A double precision floating point value giving the distance to the wall, in lattice units.

If the site is adjacent to either a wall or an iolet, there then follows data describing the distances to the wall or boundary:
  * fourteen double precision floating point values describing the distance from the site to the wall/boundary along each of the lattice directions.
  * The value is given as a fraction of the lattice vector, _not_ as a distance in lattice units.
  * If the vector does not cut the wall/boundary, then the value given is positive infinity.

### XDR Schema

Here we hope to eventually record a standardised XDR schema for the data file.
