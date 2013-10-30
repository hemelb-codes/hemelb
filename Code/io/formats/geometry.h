// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_IO_FORMATS_GEOMETRY_H
#define HEMELB_IO_FORMATS_GEOMETRY_H

#include <vector>
#include "io/formats/formats.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      /**
       * Singleton class that contains information necessary to interpret
       * a HemeLB geometry file (*.gmy).
       *
       * The file is described on the wiki http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/NewGeometryFiles
       */
      class geometry
      {
        public:
          /**
           * Get the singleton, constructing it if need be.
           * @return The single geometry instance.
           */
          static inline const geometry& Get()
          {
            if (geometry::singleton == NULL)
              geometry::singleton = new geometry();

            return *geometry::singleton;
          }

          /**
           * Magic number to identify geometry files.
           * ASCII for 'gmy', then EOF
           * Combined magic number is:
           * hex    68 6c 62 21 67 6d 79 04
           * ascii:  h  l  b  !  g  m  y EOF
           */
          enum
          {
            MagicNumber = 0x676d7904
          };

          /**
           * Version number for geometry format.
           */
          enum
          {
            VersionNumber = 4
          //!< VersionNumber
          };

          /**
           * Type codes permitted for sites
           */
          enum SiteType
          {
            SOLID = 0, //!< SOLID
            FLUID = 1
          //!< FLUID
          };

          /**
           * Type codes for the sort of iolets intersected by a link.
           */
          enum CutType
          {
            CUT_NONE = 0, //!< No intersection
            CUT_WALL = 1, //!< Intersect a wall
            CUT_INLET = 2, //!< Intersect an inlet
            CUT_OUTLET = 3
          //!< Intersect an outlet
          };

          /**
           * Type codes defining wall normal availability
           */
          enum WallNormalAvailability
          {
            WALL_NORMAL_NOT_AVAILABLE = 0, //!< WALL_NORMAL_NOT_AVAILABLE
            WALL_NORMAL_AVAILABLE = 1
          //!< WALL_NORMAL_AVAILABLE
          };

          /**
           * Number of displacements in the neighbourhood
           */
          enum
          {
            NumberOfDisplacements = 26
          };

          /**
           * The length of the preamble for the geometry file:
           *  * 1 uint for the HemeLB magic number
           *  * 1 uint for the geometry magic number
           *  * 1 uint for the version
           *  * 3 uints for the problem dimensions in blocks
           *  * 1 uint for the number of sites along one block side
           *  * 1 uint, value 0 to pad to 32 bytes
           *
           *  * 8 uints = 8 * 4  = 32
           */
          enum
          {
            PreambleLength = 32
          };

          /**
           * The length of a single header record (i.e. you have one of these
           * per block):
           *  * 1 uint for number of fluid sites
           *  * 1 uint for number of bytes occupied
           */
          enum
          {
            HeaderRecordLength = 12
          };

          /**
           * The maximum possible length of a single site's data.
           *  * 1 uint for the type
           *  * then NumberOfDisplacements:
           *    * 1 uint for the cut type
           *    * 1 uint for the inlet/outlet ID
           *    * 1 float for the cut distance
           *  * 1 uint for the wall normal availability
           *  * 3 floats for the wall normal
           */
          enum
          {
            MaxFluidSiteRecordLength = 4 + geometry::NumberOfDisplacements * (4 + 4 + 4) + 4 + 3 * 4
          };
          /**
           * Maximum length of a solid site's data. (In fact, this is THE
           * length of a solid site's data.)
           */
          enum
          {
            MaxSolidSiteRecordLength = 4
          //!< MaxSolidSiteRecordLength
          };

          /**
           * Compute the maximum possible length of a single block's data.
           * @param blockSideLength
           * @return maximum block record length in bytes
           */
          static inline unsigned int GetMaxBlockRecordLength(unsigned int blockSideLength)
          {
            return blockSideLength * blockSideLength * blockSideLength * geometry::MaxFluidSiteRecordLength;
          }

          /**
           * Compute the maximum length of a block's data, given the number of
           * fluid sites contained within it.
           * @param blockSideLength
           * @param nFluidSites
           * @return max length in bytes
           */
          static inline unsigned int GetMaxBlockRecordLength(unsigned int blockSideLength, unsigned int nFluidSites)
          {
            unsigned int nSolidSites = blockSideLength * blockSideLength * blockSideLength - nFluidSites;
            return (nFluidSites * geometry::MaxFluidSiteRecordLength + nSolidSites * geometry::MaxSolidSiteRecordLength);
          }

          /**
           * Give the displacement to a neighbouring lattice point.
           */
          typedef util::Vector3D<int> Displacement;
          typedef std::vector<Displacement> DisplacementVector;

          /**
           * Fetch the vector of displacements making up the 3D Moore
           * neighbourhood of a point, ordered as in the geometry file format.
           * @return A std::vector containing the list of displacements.
           */
          static inline const DisplacementVector& GetNeighbourhood()
          {
            return Get().displacements;
          }

        private:
          /**
           * Construct the singleton
           */
          geometry();

          /**
           * Destructor
           */
          ~geometry();

          static geometry* singleton; //! The singleton
          DisplacementVector displacements; //! The displacements of the neighbourhood

      };
    }
  }

}
#endif // HEMELB_IO_FORMATS_GEOMETRY_H
