// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_IO_FORMATS_SNAPSHOT_H
#define HEMELB_IO_FORMATS_SNAPSHOT_H

#include "io/formats/formats.h"

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      namespace snapshot
      {
        /* The format comprises a header and then a body.
         * Header is made up of (hex file position, type, description)
         * 00   uint       HemeLB magic number (see formats.h)
         * 04   uint       Snapshot magic number (see below)
         * 08   uint       Version number
         * 12   uint       Offset of the body start into the file (bytes)
         * 16   uint       Flag for simulation stability
         * 20   double     Voxel size in physical units (m)
         * 28   3 x dbl    Origin of voxel indexing in physical units (m)
         * 52   3 x uint   Bounding box minimum indices
         * 64   3 x uint   Bounding box maximum indices
         * 76   uint       Number of fluid voxels
         * Header length = 80 bytes
         *
         * Body consists of one record per voxel
         * 3 x uint     Voxel index
         * 1 x float   Pressure (mmHg)
         * 3 x float   Velocity (m/s)
         * 1 x float   von Mises stress (Pa). Note this is NaN for non-wall sites
         * Record length = 32 bytes
         *
         */

        enum
        {
          /* Identify snapshot files
           * ASCII for 'snp', then EOF
           * Combined magic number is
           * hex    68 6c 62 21 73 6e 70 04
           * ascii:  h  l  b  !  s  n  p EOF
           */
          MagicNumber = 0x736e7004
        };
        enum
        {
          VersionNumber = 2
        };
        enum
        {
          HeaderLength = 80
        };

        enum
        {
          VoxelRecordLength = 32
        };
      }
    }
  }

}
#endif // HEMELB_IO_FORMATS_SNAPSHOT_H
