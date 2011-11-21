#ifndef HEMELB_VIS_VELOCITYFIELD_H
#define HEMELB_VIS_VELOCITYFIELD_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/streaklineDrawer/NeighbouringProcessor.h"
#include "vis/streaklineDrawer/VelocitySiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class StreaklineDrawer;

      class VelocityField
      {
        public:
          VelocityField(std::map<proc_t, NeighbouringProcessor>& iNeighbouringProcessors);

          void BuildVelocityField(const geometry::LatticeData& latDat,
                                  StreaklineDrawer* iStreaklineDrawer);

          bool BlockContainsData(size_t iBlockNumber) const;

          VelocitySiteData* GetVelocitySiteData(const geometry::LatticeData& latDat,
                                                const util::Vector3D<site_t>& location);

          void GetVelocityFieldAroundPoint(const util::Vector3D<site_t> location,
                                           const geometry::LatticeData& latDat,
                                           float v[2][2][2][3]);

          void InterpolateVelocityForPoint(const float x,
                                           const float y,
                                           const float z,
                                           const float v[2][2][2][3],
                                           float interp_v[3]) const;

          void InvalidateAllCalculatedVelocities();

        private:
          // Counter to make sure the velocity field blocks are correct for the current iteration.
          site_t counter;

          VelocitySiteData& GetSiteData(site_t iBlockNumber, site_t iSiteNumber);

          void InitializeVelFieldBlock(const geometry::LatticeData& latDat, const util::Vector3D<
              site_t> location, const proc_t proc_id);

          //Vector containing VelocityFields
          std::vector<std::vector<VelocitySiteData> > velocityField;

          std::map<proc_t, NeighbouringProcessor>& neighbouringProcessors;

      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
