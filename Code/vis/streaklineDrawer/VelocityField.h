// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_STREAKLINEDRAWER_VELOCITYFIELD_H
#define HEMELB_VIS_STREAKLINEDRAWER_VELOCITYFIELD_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "lb/MacroscopicPropertyCache.h"
#include "topology/NetworkTopology.h"

#include "vis/streaklineDrawer/NeighbouringProcessor.h"
#include "vis/streaklineDrawer/VelocitySiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class VelocityField
      {
        public:
          VelocityField(std::map<proc_t, NeighbouringProcessor>& iNeighbouringProcessors
                        , const lb::MacroscopicPropertyCache& propertyCache);

          void BuildVelocityField(const geometry::LatticeData& latDat);

          bool BlockContainsData(size_t iBlockNumber) const;

          VelocitySiteData* GetVelocitySiteData(const geometry::LatticeData& latDat,
                                                const util::Vector3D<site_t>& location);

          void GetVelocityFieldAroundPoint(const util::Vector3D<site_t> location,
                                           const geometry::LatticeData& latDat,
                                           util::Vector3D<float> localVelocityField[2][2][2]);

          util::Vector3D<float>
          InterpolateVelocityForPoint(const util::Vector3D<float> position,
                                      const util::Vector3D<float> localVelocityField[2][2][2]) const;

          void InvalidateAllCalculatedVelocities();

          void UpdateLocalField(const util::Vector3D<site_t>& position, const geometry::LatticeData& latDat);

          bool NeededFromNeighbour(const util::Vector3D<site_t> location,
                                   const geometry::LatticeData& latDat,
                                   proc_t* sourceProcessor);

        private:
          void UpdateLocalField(VelocitySiteData* localVelocitySiteData, const geometry::LatticeData& latDat);

          // Counter to make sure the velocity field blocks are correct for the current iteration.
          site_t counter;

          VelocitySiteData& GetSiteData(site_t iBlockNumber, site_t iSiteNumber);

          void InitializeVelocityFieldBlock(const geometry::LatticeData& latDat,
                                            const util::Vector3D<site_t> location,
                                            const proc_t proc_id);

          // Vector containing VelocityFields
          std::vector<std::vector<VelocitySiteData> > velocityField;
          std::map<proc_t, NeighbouringProcessor>& neighbouringProcessors;
          const lb::MacroscopicPropertyCache& propertyCache;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_VELOCITYFIELD_H
