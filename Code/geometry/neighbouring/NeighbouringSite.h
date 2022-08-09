// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITE_H
#define HEMELB_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITE_H

#include "geometry/Site.h"
#include "geometry/neighbouring/NeighbouringDomain.h"

namespace hemelb
{
  namespace geometry
  {
    namespace neighbouring
    {
      class NeighbouringDomain;

//      class NeighbouringSite : public Site<NeighbouringFieldData>
//      {
//        public:
//          /**
//           * Constructor to mirror the constructor of the base class.
//           * @param localContiguousIndex
//           * @param domainData
//           */
//          NeighbouringSite(site_t localContiguousIndex, NeighbouringFieldData& domainData) :
//              Site<NeighbouringFieldData>(localContiguousIndex, domainData)
//          {
//          }
//
//          /*template<typename LatticeType>
//          inline distribn_t* GetFOld()
//          {
//            return m_fieldData->GetFOld(index * LatticeType::NUMVECTORS);
//          }
//
//          // Non-templated version of GetFOld, for when you haven't got a lattice type handy
//          inline distribn_t* GetFOld(int numvectors)
//          {
//            return m_fieldData->GetFOld(index * numvectors);
//          }*/
//
//      };

//      class ConstNeighbouringSite : public Site<const NeighbouringFieldData>
//      {
//        public:
//          /**
//           * Constructor to mirror the constructor of the base class.
//           * @param localContiguousIndex
//           * @param domainData
//           */
//          ConstNeighbouringSite(site_t localContiguousIndex,
//                                const NeighbouringFieldData& domainData) :
//              Site<const NeighbouringFieldData>(localContiguousIndex, domainData)
//          {
//          }
//      };
    }
  }
}

#endif
