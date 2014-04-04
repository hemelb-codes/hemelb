// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
