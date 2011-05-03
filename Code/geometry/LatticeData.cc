#include <map>

#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData(const bool reserveSteeringCore,
                             site_t* totalFluidSites,
                             site_t siteMins[3],
                             site_t siteMaxes[3],
                             site_t* fluidSitePerProc,
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

    void LatticeData::InitialiseNeighbourLookup(site_t** bSharedFLocationForEachProc,
                                                proc_t localRank,
                                                const unsigned int* iSiteDataForThisRank)
    {
      site_t n = -1;
      site_t lSiteIndexOnProc = 0;

      std::map<proc_t, site_t> sitesHandledPerProc;

      // Iterate over blocks in global co-ords.
      for (site_t i = 0; i < GetXSiteCount(); i += GetBlockSize())
      {
        for (site_t j = 0; j < GetYSiteCount(); j += GetBlockSize())
        {
          for (site_t k = 0; k < GetZSiteCount(); k += GetBlockSize())
          {
            n++;
            geometry::LatticeData::BlockData *map_block_p = GetBlock(n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            site_t m = -1;

            // Iterate over sites within the block.
            for (site_t site_i = i; site_i < i + GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + GetBlockSize(); site_k++)
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
                    site_t neigh_i = site_i + D3Q15::CX[l];
                    site_t neigh_j = site_j + D3Q15::CY[l];
                    site_t neigh_k = site_k + D3Q15::CZ[l];

                    if (!IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      // Set the neighbour location to the rubbish site.
                      SetNeighbourLocation(site_map, l, GetLocalFluidSiteCount()
                          * D3Q15::NUMVECTORS);
                      continue;
                    }

                    // Get the id of the processor which the neighbouring site lies on.
                    const proc_t* proc_id_p = GetProcIdFromGlobalCoords(neigh_i, neigh_j, neigh_k);

                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
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
                    else if (localRank == *proc_id_p)
                    {
                      // Pointer to the neighbour.
                      site_t contigSiteId = GetContiguousSiteId(neigh_i, neigh_j, neigh_k);

                      SetNeighbourLocation(site_map, l, contigSiteId * D3Q15::NUMVECTORS + l);

                      continue;
                    }
                    else
                    {
                      proc_t neigh_proc_index = (proc_t) *proc_id_p;

                      // This stores some coordinates.  We
                      // still need to know the site number.
                      // neigh_proc[ n ].f_data is now
                      // set as well, since this points to
                      // f_data.  Every process has data for
                      // its neighbours which say which sites
                      // on this process are shared with the
                      // neighbour.
                      site_t fluidSitesHandled = (sitesHandledPerProc.count(neigh_proc_index) > 0)
                        ? sitesHandledPerProc[neigh_proc_index]
                        : 0;

                      site_t* f_data_p =
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

    const double* LatticeData::GetNormalToWall(site_t iSiteIndex) const
    {
      return localLatDat.GetNormalToWall(iSiteIndex);
    }

    site_t LatticeData::GetXSiteCount() const
    {
      return globLatDat.GetXSiteCount();
    }
    site_t LatticeData::GetYSiteCount() const
    {
      return globLatDat.GetYSiteCount();
    }
    site_t LatticeData::GetZSiteCount() const
    {
      return globLatDat.GetZSiteCount();
    }

    site_t LatticeData::GetXBlockCount() const
    {
      return globLatDat.GetXBlockCount();
    }
    site_t LatticeData::GetYBlockCount() const
    {
      return globLatDat.GetYBlockCount();
    }
    site_t LatticeData::GetZBlockCount() const
    {
      return globLatDat.GetZBlockCount();
    }

    unsigned int LatticeData::GetLog2BlockSize() const
    {
      return globLatDat.Log2BlockSize;
    }

    site_t LatticeData::GetBlockSize() const
    {
      return globLatDat.GetBlockSize();
    }

    site_t LatticeData::GetBlockCount() const
    {
      return globLatDat.GetBlockCount();
    }

    site_t LatticeData::GetSitesPerBlockVolumeUnit() const
    {
      return globLatDat.GetSitesPerBlockVolumeUnit();
    }

    site_t LatticeData::GetBlockIdFromBlockCoords(site_t i, site_t j, site_t k) const
    {
      return globLatDat.GetBlockIdFromBlockCoords(i, j, k);
    }

    const proc_t* LatticeData::GetProcIdFromGlobalCoords(site_t siteI, site_t siteJ, site_t siteK) const
    {
      return globLatDat.GetProcIdFromGlobalCoords(siteI, siteJ, siteK);
    }

    bool LatticeData::IsValidLatticeSite(site_t i, site_t j, site_t k) const
    {
      return globLatDat.IsValidLatticeSite(i, j, k);
    }

    LatticeData::BlockData* LatticeData::GetBlock(site_t blockNumber) const
    {
      return &globLatDat.Blocks[blockNumber];
    }

    distribn_t* LatticeData::GetFOld(site_t siteNumber) const
    {
      return &localLatDat.FOld[siteNumber];
    }

    distribn_t* LatticeData::GetFNew(site_t siteNumber) const
    {
      return &localLatDat.FNew[siteNumber];
    }

    site_t LatticeData::GetLocalFluidSiteCount() const
    {
      return localLatDat.GetLocalFluidSiteCount();
    }

    LatticeData::SiteType LatticeData::GetSiteType(site_t iSiteIndex) const
    {
      return localLatDat.GetSiteType(iSiteIndex);
    }

    int LatticeData::GetBoundaryId(site_t iSiteIndex) const
    {
      return localLatDat.GetBoundaryId(iSiteIndex);
    }

    site_t LatticeData::GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const
    {
      return localLatDat.GetStreamedIndex(iSiteIndex, iDirectionIndex);
    }

    bool LatticeData::HasBoundary(site_t iSiteIndex, int iDirection) const
    {
      return localLatDat.HasBoundary(iSiteIndex, iDirection);
    }

    double LatticeData::GetCutDistance(site_t iSiteIndex, int iDirection) const
    {
      return localLatDat.GetCutDistance(iSiteIndex, iDirection);
    }

    unsigned int LatticeData::GetSiteData(site_t iSiteIndex) const
    {
      return localLatDat.mSiteData[iSiteIndex];
    }

    unsigned int LatticeData::GetContiguousSiteId(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const
    {
      return globLatDat.GetSiteData(iSiteI, iSiteJ, iSiteK);
    }

    void LatticeData::SetNeighbourLocation(site_t iSiteIndex,
                                           unsigned int iDirection,
                                           site_t iValue)
    {
      localLatDat.SetNeighbourLocation(iSiteIndex, iDirection, iValue);
    }

    void LatticeData::SetSiteCounts(site_t innerSites,
                                    site_t interCollisions[COLLISION_TYPES],
                                    site_t innerCollisions[COLLISION_TYPES],
                                    site_t sharedSites)
    {
      localLatDat.my_inner_sites = innerSites;

      for (unsigned int m = 0; m < COLLISION_TYPES; ++m)
      {
        localLatDat.my_inter_collisions[m] = interCollisions[m];
        localLatDat.my_inner_collisions[m] = innerCollisions[m];
      }

      localLatDat.SetSharedSiteCount(sharedSites);
    }

    void LatticeData::SetSiteData(site_t siteIndex, unsigned int siteData)
    {
      localLatDat.mSiteData[siteIndex] = siteData;
    }
    void LatticeData::SetWallNormal(site_t siteIndex, double normal[3])
    {
      localLatDat.SetWallNormal(siteIndex, normal);
    }
    void LatticeData::SetWallDistance(site_t siteIndex, double cutDistance[D3Q15::NUMVECTORS
        - 1])
    {
      localLatDat.SetDistanceToWall(siteIndex, cutDistance);
    }

    site_t LatticeData::GetInnerSiteCount() const
    {
      return localLatDat.my_inner_sites;
    }
    site_t LatticeData::GetInnerCollisionCount(unsigned int collisionType) const
    {
      return localLatDat.my_inner_collisions[collisionType];
    }
    site_t LatticeData::GetInterCollisionCount(unsigned int collisionType) const
    {
      return localLatDat.my_inter_collisions[collisionType];
    }

    unsigned int LatticeData::GetCollisionType(unsigned int site_data) const
    {
      return globLatDat.GetCollisionType(site_data);
    }

    void LatticeData::SwapOldAndNew()
    {
      distribn_t* temp = localLatDat.FOld;
      localLatDat.FOld = localLatDat.FNew;
      localLatDat.FNew = temp;
    }
  }
}
