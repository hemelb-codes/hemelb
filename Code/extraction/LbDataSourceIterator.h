
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
        LbDataSourceIterator(lb::MacroscopicPropertyCache& propertyCache,
                             const geometry::LatticeData& data,
                             int rank,
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
        FloatingType GetPressure() const;

        FloatingType GetDensity() const;

        /**
         * Returns the velocity at the site.
         * @return
         */
        util::Vector3D<FloatingType> GetVelocity() const;

        /**
         * Returns the shear stress at the site.
         * @return
         */
        FloatingType GetShearStress() const;

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        FloatingType GetVonMisesStress() const;

        /**
         * Returns the shear rate at the site.
         * @return shear rate
         */
        FloatingType GetShearRate() const;

        /**
         * Returns the full stress tensor at the site.
         * @return stress tensor
         */
        util::Matrix3D GetStressTensor() const;

        /**
         * Returns the traction vector at a wall site (i.e. stress tensor times surface normal).
         * @return traction vector
         */
        util::Vector3D<PhysicalStress> GetTraction() const;

        /**
         * Returns the projection of the traction vector on the tangential plane of a wall site.
         * @return projected traction vector
         */
        util::Vector3D<PhysicalStress> GetTangentialProjectionTraction() const;

        /**
         * Resets the iterator to the beginning again.
         */
        void Reset();

	/**
	 * Returns the property cache.
	 */
	lb::MacroscopicPropertyCache& GetPropertyCache();

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
        PhysicalDistance GetVoxelSize() const;

        /**
         * Returns the origin of the geometry in real, spatial units.
         *
         * @return
         */
        const PhysicalPosition& GetOrigin() const;

        /**
         * Returns true if the site at the given location is marked as a wall site
         * (i.e. one of its links intersects a wall)
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        bool IsWallSite(const util::Vector3D<site_t>& location) const;

        bool IsVesselWallSite(const util::Vector3D<site_t>& location) const;

        bool IsStentWallSite(const util::Vector3D<site_t>& location) const;

      private:
        /**
         * The cache of properties for each site, which we iterate through.
         */
        lb::MacroscopicPropertyCache& propertyCache;
        /**
         * The object containing information about the lattice.
         */
        const geometry::LatticeData& data;
        /**
         * The rank of the current process in the LB communicator.
         */
        const int rank;
        /**
         * Object capable of converting from physical to lattice units and vice versa.
         */
        const util::UnitConverter& converter;
        /**
         * Iteration variable for tracking progress through all the local fluid sites.
         */
        site_t position;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LBDATASOURCEITERATOR_H */
