// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_VIS_RAYTRACER_CLUSTERTRAVERSER_H
#define HEMELB_VIS_RAYTRACER_CLUSTERTRAVERSER_H

#include "geometry/VolumeTraverser.h"
#include "vis/rayTracer/Cluster.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      //ClusterTraverser is used to traverse the cluster
      template<typename ClusterType>
      class ClusterTraverser : public geometry::VolumeTraverser
      {
        public:
          ClusterTraverser(const ClusterType& iCluster) :
              mCluster(iCluster)
          {
          }

          virtual ~ClusterTraverser()
          {
          }

          site_t GetXCount() const
          {
            return mCluster.GetBlocksX();
          }

          site_t GetYCount() const
          {
            return mCluster.GetBlocksY();
          }

          site_t GetZCount() const
          {
            return mCluster.GetBlocksZ();
          }

        private:
          const ClusterType& mCluster;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERTRAVERSER_H
