#ifndef HEMELB_VIS_VELOCITYFIELD_H
#define HEMELB_VIS_VELOCITYFIELD_H

#include <vector>

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
          VelocityField(std::vector<NeighbouringProcessor>& iNeighbouringProcessors);

          void BuildVelocityField(const geometry::LatticeData& iLatDat,
                                  StreaklineDrawer* iStreaklineDrawer);

          bool BlockContainsData(size_t iBlockNumber);

          VelocitySiteData& GetSiteData(size_t iBlockNumber, size_t iSiteNumber);

          VelocitySiteData* velSiteDataPointer(const geometry::LatticeData& iLatDat,
                                               site_t site_i,
                                               site_t site_j,
                                               site_t site_k);

          // Counter keeps track of the number of VelSiteDatas created
          site_t counter;

          void localVelField(site_t iX,
                             site_t iY,
                             site_t iZ,
                             float v[2][2][2][3],
                             int *is_interior,
                             const geometry::LatticeData& iLatDat,
                             StreaklineDrawer* iStreaklineDrawer);

          // Private functions for initialising the velocity field.
          void initializeVelFieldBlock(const geometry::LatticeData& iLatDat,
                                       site_t site_i,
                                       site_t site_j,
                                       site_t site_k,
                                       proc_t proc_id,
                                       StreaklineDrawer* iStreaklineDrawer);

          void GetVelocityAtPoint(float x,
                                  float y,
                                  float z,
                                  float v[2][2][2][3],
                                  float interp_v[3]);

        private:

          //Vector containing VelocityFields
          std::vector<std::vector<VelocitySiteData> > mVelocityField;

          std::vector<NeighbouringProcessor>& mNeighbouringProcessors;

      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
