#ifndef HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H
#define HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H

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
         * Reads the next fluid site from the data source. Returns true if values could
         * be obtained.
         *
         * @return
         */
        bool ReadNext();

        /**
         * Returns the coordinates of the site.
         * @return
         */
        util::Vector3D<site_t> GetPosition() const;

        /**
         * Returns the pressure at the site.
         * @return
         */
        float GetPressure() const;

        /**
         * Returns the velocity at the site.
         * @return
         */
        util::Vector3D<float> GetVelocity() const;

        /**
         * Returns the shear stress at the site.
         * @return
         */
        float GetShearStress() const;

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        float GetVonMisesStress() const;

        /**
         * Returns the shear rate at the site.
         * @return shear rate
         */
        float GetShearRate() const;

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

        /**
         * Returns the origin of the geometry in real, spatial units.
         *
         * @return
         */
        const util::Vector3D<distribn_t>& GetOrigin() const;

      private:
        /**
         * The cache of properties for each site, which we iterate through.
         */
        const lb::MacroscopicPropertyCache& propertyCache;
        /**
         * The object containing information about the lattice.
         */
        const geometry::LatticeData& data;
        /**
         * Object capable of converting from physical to lattice units and vice versa.
         */
        const util::UnitConverter& converter;
        /**
         * Iteration variable for tracking progress through all the local fluid sites.
         */
        site_t position;
        /**
         * The origin of the geometry.
         */
        const util::Vector3D<distribn_t> origin;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H */
