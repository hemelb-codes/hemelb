#include <map>

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

      reader.LoadAndDecompose(&globLatDat,
                              totalFluidSites,
                              siteMins,
                              siteMaxes,
                              fluidSitePerProc,
                              bLbmParams,
                              bSimConfig,
                              lReadTime,
                              lDecomposeTime);

      int localRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &localRank);

      localLatDat.Initialise(fluidSitePerProc[localRank]);
    }

    void LatticeData::InitialiseNeighbourLookup(int ** bSharedFLocationForEachProc,
                                                int localRank,
                                                const unsigned int* iSiteDataForThisRank)
    {
      int n = -1;
      int lSiteIndexOnProc = 0;

      std::map<short int, int> sitesHandledPerProc;

      // Iterate over blocks in global co-ords.
      for (unsigned int i = 0; i < GetXSiteCount(); i += GetBlockSize())
      {
        for (unsigned int j = 0; j < GetYSiteCount(); j += GetBlockSize())
        {
          for (unsigned int k = 0; k < GetZSiteCount(); k += GetBlockSize())
          {
            n++;
            geometry::LatticeData::BlockData *map_block_p = GetBlock(n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            // Iterate over sites within the block.
            for (unsigned int site_i = i; site_i < i + GetBlockSize(); site_i++)
            {
              for (unsigned int site_j = j; site_j < j + GetBlockSize(); site_j++)
              {
                for (unsigned int site_k = k; site_k < k + GetBlockSize(); site_k++)
                {
                  // If a site is not on this process, continue.
                  m++;

                  if (localRank != map_block_p->ProcessorRankForEachBlockSite[m])
                  {
                    continue;
                  }

                  // Get site data, which is the number of the fluid site on this proc..
                  unsigned int site_map = map_block_p->site_data[m];

                  // Set neighbour location for the distribution component at the centre of
                  // this site.
                  SetNeighbourLocation(site_map, 0, site_map * D3Q15::NUMVECTORS + 0);

                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Work out positions of neighbours.
                    int neigh_i = site_i + D3Q15::CX[l];
                    int neigh_j = site_j + D3Q15::CY[l];
                    int neigh_k = site_k + D3Q15::CZ[l];

                    if (!IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // Get the id of the processor which the neighbouring site lies on.
                    int proc_id = GetProcIdFromGlobalCoords(neigh_i, neigh_j, neigh_k);

                    if (proc_id == BIG_NUMBER2 || proc_id == -1)
                    {
                      // initialize f_id to the rubbish site.
                      SetNeighbourLocation(site_map, l, GetLocalFluidSiteCount()
                          * D3Q15::NUMVECTORS);
                      continue;
                    }
                    // If on the same proc, set f_id of the
                    // current site and direction to the
                    // site and direction that it sends to.
                    // If we check convergence, the data for
                    // each site is split into that for the
                    // current and previous cycles.
                    else if (localRank == proc_id)
                    {

                      // Pointer to the neighbour.
                      unsigned int contigSiteId = GetContiguousSiteId(neigh_i, neigh_j, neigh_k);

                      SetNeighbourLocation(site_map, l, contigSiteId * D3Q15::NUMVECTORS + l);

                      continue;
                    }
                    else
                    {
                      short int neigh_proc_index = (short int) proc_id;

                      // This stores some coordinates.  We
                      // still need to know the site number.
                      // neigh_proc[ n ].f_data is now
                      // set as well, since this points to
                      // f_data.  Every process has data for
                      // its neighbours which say which sites
                      // on this process are shared with the
                      // neighbour.
                      int fluidSitesHandled =
                          (sitesHandledPerProc.count((short int) neigh_proc_index) > 0)
                            ? sitesHandledPerProc[neigh_proc_index]
                            : 0;

                      int *f_data_p =
                          &bSharedFLocationForEachProc[neigh_proc_index][fluidSitesHandled << 2];
                      f_data_p[0] = site_i;
                      f_data_p[1] = site_j;
                      f_data_p[2] = site_k;
                      f_data_p[3] = l;

                      sitesHandledPerProc[neigh_proc_index] = fluidSitesHandled + 1;
                    }
                  }

                  // This is used in Calculate BC in IO.
                  SetSiteData(site_map, iSiteDataForThisRank[lSiteIndexOnProc]);

                  if (GetCollisionType(GetSiteData(site_map)) & EDGE)
                  {
                    SetWallNormal(site_map, GetBlock(n)->wall_data[m].wall_nor);

                    SetWallDistance(site_map, GetBlock(n)->wall_data[m].cut_dist);
                  }
                  else
                  {
                    double lBigDistance[3];
                    for (unsigned int ii = 0; ii < 3; ii++)
                      lBigDistance[ii] = BIG_NUMBER;
                    SetWallNormal(site_map, lBigDistance);
                  }
                  ++lSiteIndexOnProc;
                }
              }
            }
          }
        }
      }
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
      return globLatDat.GetSitesPerBlockVolumeUnit();
    }

    unsigned int LatticeData::GetBlockIdFromBlockCoords(unsigned int i,
                                                        unsigned int j,
                                                        unsigned int k) const
    {
      return globLatDat.GetBlockIdFromBlockCoords(i, j, k);
    }

    int LatticeData::GetProcIdFromGlobalCoords(unsigned int siteI,
                                               unsigned int siteJ,
                                               unsigned int siteK) const
    {
      return globLatDat.GetProcIdFromGlobalCoords(siteI, siteJ, siteK);
    }

    bool LatticeData::IsValidLatticeSite(unsigned int i, unsigned int j, unsigned int k) const
    {
      return globLatDat.IsValidLatticeSite(i, j, k);
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

    unsigned int LatticeData::GetSiteData(unsigned int iSiteIndex) const
    {
      return localLatDat.mSiteData[iSiteIndex];
    }

    unsigned int LatticeData::GetContiguousSiteId(unsigned int iSiteI,
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

    unsigned int LatticeData::GetInnerSiteCount() const
    {
      return localLatDat.my_inner_sites;
    }
    unsigned int LatticeData::GetInnerCollisionCount(unsigned int collisionType) const
    {
      return localLatDat.my_inner_collisions[collisionType];
    }
    unsigned int LatticeData::GetInterCollisionCount(unsigned int collisionType) const
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
