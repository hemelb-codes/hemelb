#ifndef HEMELB_EXTRACTION_LBDATASOURCE_H
#define HEMELB_EXTRACTION_LBDATASOURCE_H

#include "extraction/IterableDataSource.h"
#include "geometry/LatticeData.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace extraction
  {
    class LbDataSourceIterator : public IterableDataSource
    {
      public:
        /**
         * Constructs an iterator using the LBM as a source.
         * @param propertyCache
         * @param data
         * @return
         */
        LbDataSourceIterator(const lb::MacroscopicPropertyCache& propertyCache,
                             const geometry::LatticeData& data,
                             const util::UnitConverter& converter);

        /**
         * Reads the next fluid site from the data source, obtaining its position,
         * pressure, velocity and stress. Returns true if values could be obtained.
         *
         * @param location
         * @param pressure
         * @param velocity
         * @param stress
         * @return
         */
        bool
        ReadNext(util::Vector3D<site_t>& location, float& pressure, util::Vector3D<float>& velocity, float& stress);

        /**
         * Resets the iterator to the beginning again.
         */
        void Reset();

        /**
         * Returns true iff the passed location is within the lattice.
         *
         * @param
         * @return
         */
        bool IsValidLatticeSite(const util::Vector3D<site_t>& location) const;

        /**
         * Returns true iff the given location is available on this core (i.e. if the data
         * lives on this core).
         * @return
         */
        bool IsAvailable(const util::Vector3D<site_t>& location) const;

        /**
         * Returns the real-world size of a single lattice unit.
         * @return
         */
        distribn_t GetVoxelSize() const;

      private:
        const lb::MacroscopicPropertyCache& propertyCache;
        const geometry::LatticeData& data;
        const util::UnitConverter& converter;
        site_t position;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LBDATASOURCE_H */
