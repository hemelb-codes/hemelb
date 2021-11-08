# Extracted property files

These files use XDR to encode data. A file will be created for each property requested.

We suggest the extension .xtr for these.

## Main header
The file begins with a main header (length = 60 bytes)
* uint32 - HemeLbMagicNumber
* uint32 - ExtractionMagicNumber
* uint32 - Format version number
* double - Voxel size (metres)
* double x 3 - Origin x,y,z components (metres)
* uint64 - Total number of sites
* uint32 - Field count
* uint32 - Length of the field header that follows
  
The ExtractionMagicNumber = 0x78747204
The version number is currently 4

## Field header
This header has fieldCount entries and in each one:
 * XDR string - the field name
 * uint32 - number of floats making up the field
 * double - the offset for all recorded values

## Data section
The body of the file contains a number of entries, one per timestep recorded.
Each record consists of:
 * uint64 - timestep number
 * for each output site (as many as the total number given in the main header)
  * 3x uint32 for grid position
  * for each field
    * the number of floats specified in the corresponding field header (in version 3 but not in later versions, these are doubles. In versions 4 and later, the values stored have an offset subtracted. The true value can be found from `offset + found value`.

## Offset files
These exist and will documented eventually.
