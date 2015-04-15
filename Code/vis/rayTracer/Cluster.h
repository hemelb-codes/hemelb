// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RAYTRACER_CLUSTER_H
#define HEMELB_VIS_RAYTRACER_CLUSTER_H

#include <vector>

#include "util/Vector3D.h"

#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      //The cluster structure stores data relating to the clusters
      //used by the RayTracer, in an optimal format
      //Cluster are produced by the ClusterSharedFactory
      //Caution: the data within the flow field is altered by means
      //of pointers obtained from the GetClusterSharedVoxelDataPointer
      //method
      template<typename Derived>
      class Cluster
      {
        public:
          Cluster(unsigned short xBlockCount, unsigned short yBlockCount,
                  unsigned short zBlockCount, const util::Vector3D<float>& minimalSite,
                  const util::Vector3D<float>& maximalSite,
                  const util::Vector3D<float>& minimalSiteOnMinimalBlock,
                  const util::Vector3D<site_t>& minimalBlock) :
              blocksX(xBlockCount), blocksY(yBlockCount), blocksZ(zBlockCount),
                  minSite(minimalSite), maxSite(maximalSite),
                  leastSiteOnLeastBlockInImage(minimalSiteOnMinimalBlock), minBlock(minimalBlock)
          {
          }

          unsigned int GetBlockIdFrom3DBlockLocation(
              const util::Vector3D<unsigned int>& iLocation) const
          {
            return iLocation.x * blocksY * blocksZ + iLocation.y * blocksZ + iLocation.z;
          }

          const util::Vector3D<double>* GetWallData(site_t iBlockNumber, site_t iSiteNumber) const
          {
            return ((const Derived*) (this))->DoGetWallData(iBlockNumber, iSiteNumber);
          }

          void SetWallData(site_t iBlockNumber, site_t iSiteNumber,
                           const util::Vector3D<double>& iData)
          {
            return ((Derived*) (this))->DoSetWallData(iBlockNumber, iSiteNumber, iData);
          }

          static bool NeedsWallNormals()
          {
            return Derived::DoNeedsWallNormals();
          }

          const std::vector<util::Vector3D<float> > GetCorners() const
          {
            std::vector<util::Vector3D<float> > lCorners;

            lCorners.push_back(util::Vector3D<float>(minSite.x, minSite.y, minSite.z));

            lCorners.push_back(util::Vector3D<float>(minSite.x, minSite.y, maxSite.z));

            lCorners.push_back(util::Vector3D<float>(minSite.x, maxSite.y, minSite.z));

            lCorners.push_back(util::Vector3D<float>(minSite.x, maxSite.y, maxSite.z));

            lCorners.push_back(util::Vector3D<float>(maxSite.x, minSite.y, minSite.z));

            lCorners.push_back(util::Vector3D<float>(maxSite.x, minSite.y, maxSite.z));

            lCorners.push_back(util::Vector3D<float>(maxSite.x, maxSite.y, minSite.z));

            lCorners.push_back(util::Vector3D<float>(maxSite.x, maxSite.y, maxSite.z));

            return lCorners;
          }

          /**
           * True if the cluster type requires wall normals.
           *
           * This can be overridden by deriving classes.
           * @return
           */
          static bool DoNeedsWallNormals()
          {
            return false;
          }

          unsigned short GetBlocksX() const
          {
            return blocksX;
          }

          unsigned short GetBlocksY() const
          {
            return blocksY;
          }

          unsigned short GetBlocksZ() const
          {
            return blocksZ;
          }

          const util::Vector3D<float>& GetMinSite() const
          {
            return minSite;
          }

          const util::Vector3D<float>& GetMaxSite() const
          {
            return maxSite;
          }

          const util::Vector3D<float>& GetLeastSiteOnLeastBlockInImage() const
          {
            return leastSiteOnLeastBlockInImage;
          }

          /**
           * Returns the block coordinates of the block with minimal x, y and z
           * coordinates in the cluster.
           *
           * @return The minimal block coordinates
           */
          const util::Vector3D<site_t>& GetMinBlockLocation() const
          {
            return minBlock;
          }

        private:
          /**
           * The size of the cluster in terms of the number of blocks
           */
          unsigned short blocksX;
          unsigned short blocksY;
          unsigned short blocksZ;

          /**
           * The min and maximum site location, in site units
           * relative to the centre of the lattice
           */
          util::Vector3D<float> minSite;
          util::Vector3D<float> maxSite;

          /**
           * The lowest x, y and z block location of the ClusterShared
           * in terms of site units relative to the centre location
           */
          util::Vector3D<float> leastSiteOnLeastBlockInImage;

          /**
           * The coordinates of the block with minimal x, y and z components.
           */
          util::Vector3D<site_t> minBlock;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTER_H
