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
         * Virtual destructor. We could declare this as pure virtual but that will cause a fail at
         * link-time (at least with GCC). GCC needs an object file to put the vtable in; by defining
         * this function in IterableDataSource.cc, we give it one.
         * @return
         */
        virtual ~IterableDataSource();

        /**
         * Reads the next fluid site from the data source, obtaining its position,
         * pressure, velocity and stress. Returns true if values could be obtained.
         *
         * @param position
         * @param pressure
         * @param velocity
         * @param stress
         * @return
         */
        virtual bool ReadNext(util::Vector3D<site_t>& position,
                              float& pressure,
                              util::Vector3D<float>& velocity,
                              float& stress) = 0;

        /**
         * Resets the iterator to the beginning again.
         */
        virtual void Reset() = 0;

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
        virtual distribn_t GetVoxelSize() const = 0;
    };
  }
}

#endif /* HEMELB_EXTRACTION_ITERABLEDATASOURCE_H */
