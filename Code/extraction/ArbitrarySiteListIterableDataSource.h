
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCE_H
#define HEMELB_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCE_H

#include "util/Vector3D.h"
#include "units.h"
#include "extraction/IterableDataSource.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"

namespace hemelb
{
  namespace extraction
  {
    class ArbitrarySiteListIterableDataSource : public IterableDataSource
    {
      public:
        ArbitrarySiteListIterableDataSource() {
          origin = util::Vector3D<distribn_t>(0,0,0);
        }

        /**
         * Destructor, clears up.
         * @return
         */
        ~ArbitrarySiteListIterableDataSource();

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
          //TODO: Makes this work properly!
          return origin;
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
          return origin;
        }

        /**
         * Returns true if the site at the given location is marked as a wall site
         * (i.e. one of its links intersects a wall)
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        bool IsWallSite(const util::Vector3D<site_t>& location) const {
          return false;
        }

      protected:
        std::vector<site_t> sitesWhichAreIteratedOver;
        util::Vector3D<distribn_t> origin;
    };
  }
}

#endif /* HEMELB_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCE_H */
