#ifndef HEMELB_EXTRACTION_ITERABLEDATASOURCE_H
#define HEMELB_EXTRACTION_ITERABLEDATASOURCE_H

#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace extraction
  {
    class IterableDataSource
    {
      public:
        /**
         * Reads the next fluid site from the data source, obtaining its position,
         * pressure, velocity and stress. Returns true while there are still values
         * left in the data source.
         *
         * @param position
         * @param pressure
         * @param velocity
         * @param stress
         * @return
         */
        virtual bool ReadNext(util::Vector3D<float>& position,
                              float& pressure,
                              util::Vector3D<float>& velocity,
                              float& stress) = 0;

        /**
         * Returns true iff the passed location is within the lattice.
         *
         * @param
         * @return
         */
        virtual bool IsValidLatticeSite(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns true iff the given location is available on this core (i.e. if the data
         * lives on this core).
         * @return
         */
        virtual bool IsAvailable(const util::Vector3D<site_t>& location) const = 0;

        /**
         * Returns the real-world size of a single lattice unit.
         * @return
         */
        virtual distribn_t GetVoxelSize() const;
    };
  }
}

#endif /* HEMELB_EXTRACTION_ITERABLEDATASOURCE_H */
