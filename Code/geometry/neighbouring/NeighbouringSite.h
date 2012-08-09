#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITE_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITE_H

#include "geometry/Site.h"
#include "geometry/neighbouring/NeighbouringLatticeData.h"

namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {
      class NeighbouringLatticeData;

      class NeighbouringSite : public BaseSite<NeighbouringLatticeData>
      {
        public:
          /**
           * Constructor to mirror the constructor of the base class.
           * @param localContiguousIndex
           * @param latticeData
           */
          NeighbouringSite(site_t localContiguousIndex, NeighbouringLatticeData& latticeData) :
            BaseSite<NeighbouringLatticeData> (localContiguousIndex, latticeData)
          {
          }

          template<typename LatticeType>
          inline distribn_t* GetFOld()
          {
            return latticeData.GetFOld(index * LatticeType::NUMVECTORS);
          }

          // Non-templated version of GetFOld, for when you haven't got a lattice type handy
          inline distribn_t* GetFOld(int numvectors)
          {
            return latticeData.GetFOld(index * numvectors);
          }

      };

      class ConstNeighbouringSite : public BaseSite<const NeighbouringLatticeData>
      {
        public:
          /**
           * Constructor to mirror the constructor of the base class.
           * @param localContiguousIndex
           * @param latticeData
           */
          ConstNeighbouringSite(site_t localContiguousIndex, const NeighbouringLatticeData& latticeData) :
            BaseSite<const NeighbouringLatticeData> (localContiguousIndex, latticeData)
          {
          }
      };
    }
  }
}

#endif
