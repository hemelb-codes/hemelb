/*
 * ArbitrarySiteListIterableDataSource.h
 *
 *  Created on: 26 Jun 2012
 *      Author: derek
 */

#ifndef HEMELB_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCE_H
#define HEMELB_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCE_H

#include "util/Vector3D.h"
#include "units.h"
#include "extraction/IterableDataSource.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "geometry/LatticeData.h"
#include "geometry/Site.h"

namespace hemelb
{
  namespace extraction
  {
    class ArbitrarySiteListIterableDataSource : public IterableDataSource
    {
      public:
        /**
         * Virtual destructor. We could declare this as pure virtual but that will cause a fail at
         * link-time (at least with GCC). GCC needs an object file to put the vtable in; by defining
         * this function in IterableDataSource.cc, we give it one.
         * @return
         */
        ~ArbitrarySiteListIterableDataSource();

        void SetManager(geometry::neighbouring::NeighbouringDataManager *manager)
        {
          /* NeighbouringLatticeData is an instance of BaseSite. Do we want to store sites or indices here?
           * Doing full sites for now, may speed up computation as the iteration pattern is more
           * predictable? */
          iterableData = manager;
        }

        /**
         * Returns the velocity at the site.
         * @return
         */
        distribn_t GetVelocityRelativeToNormal(geometry::neighbouring::NeighbouringDataManager *NDM, distribn_t *normal)
        {
//          distribn_t velocity[3];
//          distribn_t total_velocity[3];
//          distribn_t density;
          return 1.0;

//          /* Apply CalcDensityAndVelocity to extract velocities and add them all up.
//           * We're not (yet) using weights or normalisation here. */
//          for (int i = 0; i < iterableData.size(); i++)
//          {
//            std::vector<distribn_t> f_data = NLD->GetDistribution(iterableData[i].GetIndex());
//            LatticeType::CalculateDensityAndVelocity(f_data, &density, &velocity[0], &velocity[1], &velocity[2]);
//            total_velocity[0] += velocity[0];
//            total_velocity[1] += velocity[1];
//            total_velocity[2] += velocity[2];
//          }
//
//          /* Dot product of the total velocity with the boundary normal. */
//          return (velocity[0] * normal[0]) + (velocity[1] * normal[1]) + (velocity[2] * normal[2]);
        }

        /**
         * Reads the next fluid site from the data source. Its position,
         * pressure, velocity, etc are available from other functions.
         * Returns true if values could be obtained.
         *
         * @return
         */
        bool ReadNext() {
          return false;
        }

        /**
         * Returns the coordinates of the site.
         * @return
         */
        util::Vector3D<site_t> GetPosition() const {
          return util::Vector3D<site_t>(0, 0, 0);
        }

        /**
         * Returns the pressure at the site.
         * @return
         */
        float GetPressure() const {
          return 0.0;
        }

        /**
         * Returns the velocity at the site.
         * @return
         */
        util::Vector3D<float> GetVelocity() const {
          return util::Vector3D<site_t>(0.0, 0.0, 0.0);
        }

        /**
         * Returns the shear stress at the site.
         * @return
         */
        float GetShearStress() const {
          return 0.0;
        }

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        float GetVonMisesStress() const {
          return 0.0;
        }

        /**
         * Returns the shear rate at the site.
         * @return
         */
        float GetShearRate() const {
          return 0.0;
        }

        /**
         * Resets the iterator to the beginning again.
         */
        void Reset() {
          return;
        }

        /**
         * Returns true iff the passed location is within the lattice.
         *
         * @param
         * @return
         */
        bool IsValidLatticeSite(const util::Vector3D<site_t>& location) const {
          return false;
        }

        /**
         * Returns true iff the given location is available on this core (i.e. if the data
         * lives on this core).
         * @return
         */
        bool IsAvailable(const util::Vector3D<site_t>& location) const {
          return false;
        }

        /**
         * Returns the real-world size of a single lattice unit.
         * @return
         */
        distribn_t GetVoxelSize() const {
          return 0.0;
        }

        /**
         * Returns the origin of the geometry in real units.
         * @return
         */
        const util::Vector3D<distribn_t>& GetOrigin() const {
          return util::Vector3D<site_t>(0, 0, 0);
        }

        /**
         * Returns true if the site at the given location is marked as an edge site
         * (i.e. one of its links intersects a wall)
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        bool IsEdgeSite(const util::Vector3D<site_t>& location) const {
          return false;
        }

      protected:
        geometry::neighbouring::NeighbouringDataManager *iterableData;
    };
  }
}

#endif /* HEMELB_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCE_H */
