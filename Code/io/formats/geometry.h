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
           * Magic number to identify snapshot files.
           * ASCII for 'gmy', then EOF
           * Combined magic number is:
           * hex    68 6c 62 21 67 6d 79 04
           * ascii:  h  l  b  !  s  n  p EOF
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
            VersionNumber = 2//!< VersionNumber
          };

          /**
           * The length of the preamble for the geometry file:
           *  * 1 uint for the HemeLB magic number
           *  * 1 uint for the geometry magic number
           *  * 1 uint for the version
           *  * 3 uints for the problem dimensions in blocks
           *  * 1 uint for the number of sites along one block side
           *  * 1 double for the voxel size
           *  * 3 doubles for the coords of the origin
           *  * 1 uint, value 0 to pad to 64 bytes
           *
           *  * 8 uints, 4 doubles = 8 * 4 + 4 * 8 = 64
           */
          static const unsigned PreambleLength = 64;

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
