
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

      class NeighbouringSite : public Site<NeighbouringLatticeData>
      {
        public:
          /**
           * Constructor to mirror the constructor of the base class.
           * @param localContiguousIndex
           * @param latticeData
           */
          NeighbouringSite(site_t localContiguousIndex, NeighbouringLatticeData& latticeData) :
              Site<NeighbouringLatticeData>(localContiguousIndex, latticeData)
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

      class ConstNeighbouringSite : public Site<const NeighbouringLatticeData>
      {
        public:
          /**
           * Constructor to mirror the constructor of the base class.
           * @param localContiguousIndex
           * @param latticeData
           */
          ConstNeighbouringSite(site_t localContiguousIndex, const NeighbouringLatticeData& latticeData) :
              Site<const NeighbouringLatticeData>(localContiguousIndex, latticeData)
          {
          }
      };
    }
  }
}

#endif
