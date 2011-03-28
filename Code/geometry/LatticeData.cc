#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData(const bool reserveSteeringCore,
                             int* totalFluidSites,
                             unsigned int siteMins[3],
                             unsigned int siteMaxes[3],
                             unsigned int* fluidSitePerProc,
                             lb::LbmParameters* bLbmParams,
                             SimConfig* bSimConfig,
                             double* lReadTime,
                             double* lDecomposeTime) :
      localLatDat(), globLatDat()

    {
      GeometryReader reader(reserveSteeringCore);

      reader.LoadAndDecompose(&globLatDat, totalFluidSites, siteMins, siteMaxes, fluidSitePerProc,
                              bLbmParams, bSimConfig, lReadTime, lDecomposeTime);

      int localRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &localRank);

      localLatDat.Initialise(fluidSitePerProc[localRank]);
    }

    const double* LatticeData::GetNormalToWall(int iSiteIndex) const
    {
      return localLatDat.GetNormalToWall(iSiteIndex);
    }

    unsigned int LatticeData::GetXSiteCount() const
    {
      return globLatDat.GetXSiteCount();
    }
    unsigned int LatticeData::GetYSiteCount() const
    {
      return globLatDat.GetYSiteCount();
    }
    unsigned int LatticeData::GetZSiteCount() const
    {
      return globLatDat.GetZSiteCount();
    }

    unsigned int LatticeData::GetXBlockCount() const
    {
      return globLatDat.GetXBlockCount();
    }
    unsigned int LatticeData::GetYBlockCount() const
    {
      return globLatDat.GetYBlockCount();
    }
    unsigned int LatticeData::GetZBlockCount() const
    {
      return globLatDat.GetZBlockCount();
    }

    unsigned int LatticeData::GetLog2BlockSize() const
    {
      return globLatDat.Log2BlockSize;
    }

    unsigned int LatticeData::GetBlockSize() const
    {
      return globLatDat.GetBlockSize();
    }

    unsigned int LatticeData::GetBlockCount() const
    {
      return globLatDat.GetBlockCount();
    }

    unsigned int LatticeData::GetSitesPerBlockVolumeUnit() const
    {
      return globLatDat.SitesPerBlockVolumeUnit;
    }

    unsigned int LatticeData::GetBlockIdFromBlockCoords(unsigned int i,
                                                        unsigned int j,
                                                        unsigned int k) const
    {
      return globLatDat.GetBlockIdFromBlockCoords(i, j, k);
    }

    int* LatticeData::GetProcIdFromGlobalCoords(unsigned int siteI,
                                                unsigned int siteJ,
                                                unsigned int siteK) const
    {
      return globLatDat.GetProcIdFromGlobalCoords(siteI, siteJ, siteK);
    }

    LatticeData::BlockData* LatticeData::GetBlock(unsigned int blockNumber) const
    {
      return &globLatDat.Blocks[blockNumber];
    }

    double* LatticeData::GetFOld(unsigned int siteNumber) const
    {
      return &localLatDat.FOld[siteNumber];
    }

    double* LatticeData::GetFNew(unsigned int siteNumber) const
    {
      return &localLatDat.FNew[siteNumber];
    }

    unsigned int LatticeData::GetLocalFluidSiteCount() const
    {
      return localLatDat.GetLocalFluidSiteCount();
    }

    LatticeData::SiteType LatticeData::GetSiteType(int iSiteIndex) const
    {
      return localLatDat.GetSiteType(iSiteIndex);
    }

    int LatticeData::GetBoundaryId(int iSiteIndex) const
    {
      return localLatDat.GetBoundaryId(iSiteIndex);
    }

    int LatticeData::GetStreamedIndex(int iSiteIndex, int iDirectionIndex) const
    {
      return localLatDat.GetStreamedIndex(iSiteIndex, iDirectionIndex);
    }

    bool LatticeData::HasBoundary(int iSiteIndex, int iDirection) const
    {
      return localLatDat.HasBoundary(iSiteIndex, iDirection);
    }

    double LatticeData::GetCutDistance(int iSiteIndex, int iDirection) const
    {
      return localLatDat.GetCutDistance(iSiteIndex, iDirection);
    }

    unsigned int* LatticeData::GetSiteData(unsigned int iSiteIndex) const
    {
      return &localLatDat.mSiteData[iSiteIndex];
    }

    const unsigned int* LatticeData::GetSiteData(unsigned int iSiteI,
                                                 unsigned int iSiteJ,
                                                 unsigned int iSiteK) const
    {
      return globLatDat.GetSiteData(iSiteI, iSiteJ, iSiteK);
    }

    void LatticeData::SetNeighbourLocation(unsigned int iSiteIndex,
                                           unsigned int iDirection,
                                           unsigned int iValue)
    {
      localLatDat.SetNeighbourLocation(iSiteIndex, iDirection, iValue);
    }

    void LatticeData::SetSiteCounts(unsigned int innerSites,
                                    unsigned int interCollisions[COLLISION_TYPES],
                                    unsigned int innerCollisions[COLLISION_TYPES],
                                    unsigned int sharedSites)
    {
      localLatDat.my_inner_sites = innerSites;

      for (unsigned int m = 0; m < COLLISION_TYPES; ++m)
      {
        localLatDat.my_inter_collisions[m] = interCollisions[m];
        localLatDat.my_inner_collisions[m] = innerCollisions[m];
      }

      localLatDat.SetSharedSiteCount(sharedSites);
    }

    void LatticeData::SetSiteData(unsigned int siteIndex, unsigned int siteData)
    {
      localLatDat.mSiteData[siteIndex] = siteData;
    }
    void LatticeData::SetWallNormal(unsigned int siteIndex, double normal[3])
    {
      localLatDat.SetWallNormal(siteIndex, normal);
    }
    void LatticeData::SetWallDistance(unsigned int siteIndex, double cutDistance[D3Q15::NUMVECTORS
        - 1])
    {
      localLatDat.SetDistanceToWall(siteIndex, cutDistance);
    }

    unsigned int LatticeData::GetInnerSiteCount()
    {
      return localLatDat.my_inner_sites;
    }
    unsigned int LatticeData::GetInnerCollisionCount(unsigned int collisionType)
    {
      return localLatDat.my_inner_collisions[collisionType];
    }
    unsigned int LatticeData::GetInterCollisionCount(unsigned int collisionType)
    {
      return localLatDat.my_inter_collisions[collisionType];
    }

    unsigned int LatticeData::GetCollisionType(unsigned int site_data) const
    {
      return globLatDat.GetCollisionType(site_data);
    }

    void LatticeData::SwapOldAndNew()
    {
      double *temp = localLatDat.FOld;
      localLatDat.FOld = localLatDat.FNew;
      localLatDat.FNew = temp;
    }
  }
}
