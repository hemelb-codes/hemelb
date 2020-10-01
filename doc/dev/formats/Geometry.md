# Geometry File Format

This page describes the file format used to describe HemeLB lattice geometries.

The files should have the extension .gmy

## Description

This file describes the problem domain as a series of blocks, each block being an identically sized cubic subsection of the problem domain.

The file is binary data using the [XDR standard](http://tools.ietf.org/html/rfc4506). This is read/written using standard libraries, in both C, using [<rpc/xdr.h>](http://linux.die.net/man/3/xdr), and in Python, using [xdrlib](http://docs.python.org/library/xdrlib.html).

This file describes the problem domain as a series of  blocks, each block being an identically sized cubic subsection of the problem domain.  The data for each block is compressed with zlib (preamble and headers are uncompressed).


### Preamble

The file begins with an uncompressed preamble. This consists of:

* An unsigned int giving the HemeLB magic number (0x686c6221; ASCII for 'hlb!'; see source:Code/io/formats/formats.h)
* An unsigned int giving the Geometry magic number (0x676d7904; ASCII for 'gmy' then EOF; see source:Code/io/formats/geometry.h)
* An unsigned int giving the version number
* Three unsigned integers representing the x, y, and z size of the problem domain, in blocks.
* An unsigned integer representing the number of lattice sites along the side of a block. (So the total sites in a block is the cube of this.)
* An unsigned int with the value zero to pad the length of the preamble to 32 bytes

### Block Headers
The header is uncompressed and consists of an array of triples of unsigned integers describing each block. 

Blocks are striped with the z coordinate changing most frequently (i.e. the z coordinate represents the least significant part of the block index) then y and x slowest.

Each triple is:

* An unsigned integer giving the number of fluid lattice sites in the block.
  * This does not have to be the same as the cube of the side length of a block in lattice sites, because some of the sites in the block might be "solid" sites which are not simulated in the Lattice Boltzmann.
* An unsigned integer giving the total length, in bytes, of the _compressed_ data representing a given block
  * This is useful to allow one to jump to the data for a given block (by computing the partial sum of these offsets up to the block index).
* An unsigned integer giving the total length, in bytes, of the _uncompressed_ data representing a given block
### Block data

* Only if a block contains at least one fluid site is data given for the block. 
* If a block has no fluid sites, both the number of sites in the block, and the length of the data for the block, as given in the block headers array, will be zero. 
* A parser must use this knowledge so that the index of the block for a given sequence of block data is established correctly.

* The data for a block consists of data for each site in the block. 
* Sites are striped with the z coordinate changing most frequently, i.e. the z coordinate represents the least significant part of the site index within the block. 
* For each site, the data given, and the length of the data given, depends on the type of the site. 
* The data is given site by site, one after the other, without padding. Because of this, it is not possible to find the location in the block data of the data for a particular site without reading all the preceding sites within the block
* This data is compressed with [zlib](http://zlib.net/manual.html) on a per-block basis.

### Sites
For each site the following is given, in order:

* An unsigned integer indicating solid (0) or fluid (1).
* If solid, this is all
* If fluid, there follows for every link to sites in the 3D Moore neighbourhood (i.e. for all 26 sites).
* The ordering of links is as given in source:Code/io/formats/geometry.cc
* The data given per link is:
 * An unsigned integer code indicating the type of boundary crossed
 * Further information may follow each type, depending on the type
  * 0 - NONE
   * no data
  * 1 - WALL
   * a single precision floating point number giving the distance to the boundary, as a fraction of the lattice vector, _not_ as a distance in lattice units.
  * 2 - INLET
   * an unsigned integer giving the an index to the inlet array specified in the [XML Config File](XMLConfigFile/)
   * a single precision floating point number giving the distance to the boundary, as a fraction of the lattice vector, _not_ as a distance in lattice units.
  * 3 - OUTLET
   * an unsigned integer giving the an index to the outlet array specified in the [XML Config File](XMLConfigFile/)
   * a single precision floating point number giving the distance to the boundary, as a fraction of the lattice vector, _not_ as a distance in lattice units.
* An unsigned integer indicating whether a wall normal is available at the current fluid site (yes=1, no=0).
 * A wall normal must be available if one or more links intersect a wall (see previous site links description). 
* If normal available, three single precision floating point numbers specifying x, y, and z component of the wall normal.


## XDR Schema

Here we hope to eventually record a standardised XDR schema for the data file.
