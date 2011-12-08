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
          Cluster(unsigned short xBlockCount,
                  unsigned short yBlockCount,
                  unsigned short zBlockCount,
                  const util::Vector3D<float>& minimalSite,
                  const util::Vector3D<float>& maximalSite,
                  const util::Vector3D<float>& minimalSiteOnMinimalBlock) :
            blocksX(xBlockCount), blocksY(yBlockCount), blocksZ(zBlockCount), minSite(minimalSite),
                maxSite(maximalSite), leastSiteOnLeastBlockInImage(minimalSiteOnMinimalBlock)

          {
            SiteData.resize(blocksX * blocksY * blocksZ);
          }

          unsigned int GetBlockIdFrom3DBlockLocation(const util::Vector3D<unsigned int>& iLocation) const
          {
            return iLocation.x * blocksY * blocksZ + iLocation.y * blocksZ + iLocation.z;
          }

          void ResizeVectorsForBlock(site_t iBlockNumber, site_t numberOfSites)
          {
            SiteData[iBlockNumber].resize(numberOfSites, SiteData_t(-1.0F));
            ((Derived*) (this))->DoResizeVectorsForBlock(iBlockNumber, numberOfSites);
          }

          //Returns true if there is site data for a given block
          bool BlockContainsSites(site_t iBlockNumber) const
          {
            return !SiteData[iBlockNumber].empty();
          }

          const SiteData_t& GetSiteData(site_t iBlockNumber, site_t iSiteNumber) const
          {
            return SiteData[iBlockNumber][iSiteNumber];
          }

          SiteData_t& GetSiteData(site_t iBlockNumber, site_t iSiteNumber)
          {
            return SiteData[iBlockNumber][iSiteNumber];
          }

          const double* GetWallData(site_t iBlockNumber, site_t iSiteNumber) const
          {
            return ((const Derived*) (this))->DoGetWallData(iBlockNumber, iSiteNumber);
          }

          void SetWallData(site_t iBlockNumber, site_t iSiteNumber, const double* const iData)
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

        private:
          /**
           * Default version of this function doesn't need to do anything extra
           *
           * @param iBlockNumber
           * @param iSize
           */
          void DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
          {

          }

          //Stores the size of the cluster in terms of the number of blocks
          unsigned short blocksX;
          unsigned short blocksY;
          unsigned short blocksZ;

          //The min and maximum site location, in site units
          //relative to the centre of the lattice
          util::Vector3D<float> minSite;
          util::Vector3D<float> maxSite;

          //Stores the lowest x, y and z block location of the ClusterShared
          //in terms of site units relative to the centre location
          util::Vector3D<float> leastSiteOnLeastBlockInImage;

          //    public:
          std::vector<std::vector<SiteData_t> > SiteData;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTER_H
