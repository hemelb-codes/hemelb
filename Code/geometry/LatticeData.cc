#include <map>
#include <limits>

#include "log/Logger.h"
#include "topology/NetworkTopology.h"
#include "geometry/LatticeData.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData()
    {
      neighbourIndexFromProcRank = NULL;
      f_recv_iv = NULL;
    }

    LatticeData::~LatticeData()
    {
      delete[] neighbourIndexFromProcRank;
      delete[] f_recv_iv;
    }

    LatticeData* LatticeData::Load(const bool reserveSteeringCore,
                                   std::string& dataFilePath,
                                   lb::LbmParameters* bLbmParams,
                                   reporting::Timers &timings)
    {
      // Use a reader to read in the file.
      GeometryReader reader(reserveSteeringCore);

      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Loading file and decomposing geometry.");
      GlobalLatticeData* globalLattice = reader.LoadAndDecompose(bLbmParams, dataFilePath, timings);

      globalLattice->CollectFluidSiteDistribution();

      // Count the fluid sites on the local processor.
      proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();
      proc_t size = topology::NetworkTopology::Instance()->GetProcessorCount();

      //TODO this is a total hack just for now.
      site_t localMins[3];
      site_t localMaxes[3];
      localMins[0] = std::numeric_limits<site_t>::max();
      localMins[1] = std::numeric_limits<site_t>::max();
      localMins[2] = std::numeric_limits<site_t>::max();
      localMaxes[0] = 0;
      localMaxes[1] = 0;
      localMaxes[2] = 0;

      for (site_t siteI = 0; siteI < globalLattice->GetXSiteCount(); ++siteI)
      {
        for (site_t siteJ = 0; siteJ < globalLattice->GetYSiteCount(); ++siteJ)
        {
          for (site_t siteK = 0; siteK < globalLattice->GetZSiteCount(); ++siteK)
          {
            const proc_t
                * procId = globalLattice->GetProcIdFromGlobalCoords(util::Vector3D<site_t>(siteI,
                                                                                           siteJ,
                                                                                           siteK));
            if (procId == NULL || *procId != localRank)
            {
              continue;
            }
            else
            {
              localMins[0] = hemelb::util::NumericalFunctions::min(localMins[0], siteI);
              localMins[1] = hemelb::util::NumericalFunctions::min(localMins[1], siteJ);
              localMins[2] = hemelb::util::NumericalFunctions::min(localMins[2], siteK);
              localMaxes[0] = hemelb::util::NumericalFunctions::max(localMaxes[0], siteI);
              localMaxes[1] = hemelb::util::NumericalFunctions::max(localMaxes[1], siteJ);
              localMaxes[2] = hemelb::util::NumericalFunctions::max(localMaxes[2], siteK);
            }
          }
        }
      }

      site_t siteMins[3], siteMaxes[3];
      MPI_Allreduce(localMins, siteMins, 3, MpiDataType<site_t> (), MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, siteMaxes, 3, MpiDataType<site_t> (), MPI_MAX, MPI_COMM_WORLD);

      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising localLatDat.");
      LocalLatticeData
          * localLattice =
              new LocalLatticeData(globalLattice->GetFluidSiteCountOnProc(topology::NetworkTopology::Instance()->GetLocalRank()));

      for (unsigned direction = 0; direction < 3; direction++)
      {
        globalLattice->globalSiteMins[direction] = siteMins[direction];
        globalLattice->globalSiteMaxes[direction] = siteMaxes[direction];
      }

      LatticeData* lattice = new LatticeData(localLattice, globalLattice);

      lattice->totalFluidSites = 0;
      for (proc_t ii = 0; ii < (proc_t) size; ++ii)
      {
        lattice->totalFluidSites += globalLattice->fluidSitesOnEachProcessor[ii];
      }

      return lattice;
    }

    LatticeData::LatticeData(LocalLatticeData* localLattice, GlobalLatticeData* globalLattice) :
      localLatDat(*localLattice), globLatDat(*globalLattice)
    {
      Initialise();
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
            geometry::BlockData *map_block_p = GetBlock(n);

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
                  localLatDat.SetNeighbourLocation(site_map, 0, site_map * D3Q15::NUMVECTORS + 0);

                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Work out positions of neighbours.
                    site_t neighbourI = site_i + D3Q15::CX[l];
                    site_t neighbourJ = site_j + D3Q15::CY[l];
                    site_t neighbourK = site_k + D3Q15::CZ[l];

                    if (!IsValidLatticeSite(neighbourI, neighbourJ, neighbourK))
                    {
                      // Set the neighbour location to the rubbish site.
                      localLatDat.SetNeighbourLocation(site_map,
                                                       l,
                                                       GetLocalFluidSiteCount() * D3Q15::NUMVECTORS);
                      continue;
                    }

                    // Get the id of the processor which the neighbouring site lies on.
                    const proc_t * proc_id_p =
                        GetProcIdFromGlobalCoords(util::Vector3D<site_t>(neighbourI,
                                                                         neighbourJ,
                                                                         neighbourK));

                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                    {
                      // initialize f_id to the rubbish site.
                      localLatDat.SetNeighbourLocation(site_map,
                                                       l,
                                                       GetLocalFluidSiteCount() * D3Q15::NUMVECTORS);
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
                      site_t contigSiteId = GetContiguousSiteId(neighbourI, neighbourJ, neighbourK);

                      localLatDat.SetNeighbourLocation(site_map,
                                                       l,
                                                       contigSiteId * D3Q15::NUMVECTORS + l);

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

                  if (GetCollisionType(GetSiteData(site_map)) & (EDGE|INLET|OUTLET))
                  {
                    SetWallNormal(site_map, GetBlock(n)->wall_data[m].wall_nor);
                    SetWallDistance(site_map, GetBlock(n)->wall_data[m].cut_dist);
                  }
                  else
                  {
                    double lBigDistance[14];
                    for (unsigned int ii = 0; ii < 14; ii++)
                      lBigDistance[ii] = NO_VALUE;
                    SetWallDistance(site_map, lBigDistance);
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
    distribn_t LatticeData::GetVoxelSize() const
    {
      return globLatDat.GetVoxelSize();
    }
    distribn_t LatticeData::GetXOrigin() const
    {
      return globLatDat.GetXOrigin();
    }
    distribn_t LatticeData::GetYOrigin() const
    {
      return globLatDat.GetYOrigin();
    }
    distribn_t LatticeData::GetZOrigin() const
    {
      return globLatDat.GetZOrigin();
    }

    unsigned int LatticeData::GetLog2BlockSize() const
    {
      return globLatDat.log2BlockSize;
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

    const proc_t* LatticeData::GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const
    {
      return globLatDat.GetProcIdFromGlobalCoords(globalSiteCoords);
    }

    bool LatticeData::IsValidBlock(site_t i, site_t j, site_t k) const
    {
      return globLatDat.IsValidBlock(i, j, k);
    }

    bool LatticeData::IsValidLatticeSite(site_t i, site_t j, site_t k) const
    {
      return globLatDat.IsValidLatticeSite(i, j, k);
    }

    BlockData* LatticeData::GetBlock(site_t blockNumber) const
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

    bool LatticeData::HasBoundary(const site_t iSiteIndex, const int iDirection) const
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

    unsigned int LatticeData::GetContiguousSiteId(site_t siteI, site_t siteJ, site_t siteK) const
    {
      return globLatDat.GetSiteData(siteI, siteJ, siteK);
    }

    const util::Vector3D<site_t> LatticeData::GetGlobalCoords(site_t blockNumber,
                                                              const util::Vector3D<site_t>& localSiteCoords) const
    {
      return globLatDat.GetGlobalCoords(blockNumber, localSiteCoords);
    }

    void LatticeData::SetSiteData(site_t siteIndex, unsigned int siteData)
    {
      localLatDat.mSiteData[siteIndex] = siteData;
    }
    void LatticeData::SetWallNormal(site_t siteIndex, double normal[3])
    {
      localLatDat.SetWallNormal(siteIndex, normal);
    }
    void LatticeData::SetWallDistance(site_t siteIndex, double cutDistance[D3Q15::NUMVECTORS - 1])
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

    void LatticeData::CountCollisionTypes(const unsigned int * lThisRankSiteData)
    {
      site_t innerSites = 0;

      site_t interCollisions[COLLISION_TYPES];
      site_t innerCollisions[COLLISION_TYPES];

      for (unsigned int m = 0; m < COLLISION_TYPES; m++)
      {
        interCollisions[m] = 0;
        innerCollisions[m] = 0;
      }

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      totalSharedFs = 0;

      int lSiteIndexOnProc = 0;

      site_t n = -1;

      // Iterate over all blocks in site units
      for (site_t i = 0; i < GetXSiteCount(); i += GetBlockSize())
      {
        for (site_t j = 0; j < GetYSiteCount(); j += GetBlockSize())
        {
          for (site_t k = 0; k < GetZSiteCount(); k += GetBlockSize())
          {
            BlockData * map_block_p = GetBlock(++n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            site_t m = -1;

            // Iterate over all sites within the current block.
            for (site_t site_i = i; site_i < i + GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + GetBlockSize(); site_k++)
                {
                  m++;
                  // If the site is not on this processor, continue.
                  if (netTop->GetLocalRank() != map_block_p->ProcessorRankForEachBlockSite[m])
                  {
                    continue;
                  }

                  bool lIsInnerSite = true;

                  // Iterate over all direction vectors.
                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Find the neighbour site co-ords in this direction.
                    site_t neigh_i = site_i + D3Q15::CX[l];
                    site_t neigh_j = site_j + D3Q15::CY[l];
                    site_t neigh_k = site_k + D3Q15::CZ[l];

                    if (!IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // Find the processor Id for that neighbour.
                    const proc_t
                        * proc_id_p = GetProcIdFromGlobalCoords(util::Vector3D<site_t>(neigh_i,
                                                                                       neigh_j,
                                                                                       neigh_k));

                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                    // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                    // in lbmReadConfig in io.cc.
                    if (proc_id_p == NULL || netTop->GetLocalRank() == *proc_id_p || *proc_id_p
                        == (BIG_NUMBER2))
                    {
                      continue;
                    }

                    lIsInnerSite = false;

                    // The first time, net_neigh_procs = 0, so
                    // the loop is not executed.
                    bool flag = true;

                    // Iterate over neighbouring processors until we find the one with the
                    // neighbouring site on it.
                    proc_t lNeighbouringProcs = (proc_t) neighbouringProcs.size();
                    for (proc_t mm = 0; mm < lNeighbouringProcs && flag; mm++)
                    {
                      // Check whether the rank for a particular neighbour has already been
                      // used for this processor.  If it has, set flag to zero.
                      NeighbouringProcessor* neigh_proc_p = &neighbouringProcs[mm];

                      // If ProcessorRankForEachBlockSite is equal to a neigh_proc that has alredy been listed.
                      if (*proc_id_p == neigh_proc_p->Rank)
                      {
                        flag = false;
                        ++neigh_proc_p->SharedFCount;
                        ++totalSharedFs;
                      }
                    }
                    // If flag is 1, we didn't find a neighbour-proc with the neighbour-site on it
                    // so we need a new neighbouring processor.
                    if (flag)
                    {
                      // Store rank of neighbour in >neigh_proc[neigh_procs]
                      NeighbouringProcessor lNewNeighbour;
                      lNewNeighbour.SharedFCount = 1;
                      lNewNeighbour.Rank = *proc_id_p;
                      neighbouringProcs.push_back(lNewNeighbour);
                      ++totalSharedFs;
                    }
                  }

                  // Set the collision type data. map_block site data is renumbered according to
                  // fluid site numbers within a particular collision type.

                  int l = -1;

                  switch (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc]))
                  {
                    case FLUID:
                      l = 0;
                      break;
                    case EDGE:
                      l = 1;
                      break;
                    case INLET:
                      l = 2;
                      break;
                    case OUTLET:
                      l = 3;
                      break;
                    case (INLET | EDGE):
                      l = 4;
                      break;
                    case (OUTLET | EDGE):
                      l = 5;
                      break;
                  }

                  ++lSiteIndexOnProc;

                  if (lIsInnerSite)
                  {
                    ++innerSites;

                    if (l == 0)
                    {
                      map_block_p->site_data[m] = (unsigned int) innerCollisions[l];
                    }

                    ++innerCollisions[l];
                  }
                  else
                  {
                    ++interCollisions[l];
                  }
                }
              }
            }
          }
        }
      }

      site_t collision_offset[2][COLLISION_TYPES];
      // Calculate the number of each type of collision.
      collision_offset[0][0] = 0;

      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        collision_offset[0][l] = collision_offset[0][l - 1] + innerCollisions[l - 1];
      }
      collision_offset[1][0] = innerSites;
      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        collision_offset[1][l] = collision_offset[1][l - 1] + interCollisions[l - 1];
      }

      lSiteIndexOnProc = 0;

      site_t innerColsPassed[COLLISION_TYPES];
      site_t interColsPassed[COLLISION_TYPES];

      for (unsigned int ii = 0; ii < COLLISION_TYPES; ++ii)
      {
        innerColsPassed[ii] = 0;
        interColsPassed[ii] = 0;
      }

      // Iterate over blocks
      n = -1;

      // Iterate over all blocks in site units
      for (site_t i = 0; i < GetXSiteCount(); i += GetBlockSize())
      {
        for (site_t j = 0; j < GetYSiteCount(); j += GetBlockSize())
        {
          for (site_t k = 0; k < GetZSiteCount(); k += GetBlockSize())
          {
            BlockData * map_block_p = GetBlock(++n);

            // If we are in a block of solids, continue.
            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            // Iterate over sites within the block.
            site_t m = -1;

            // Iterate over all sites within the current block.
            for (site_t site_i = i; site_i < i + GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + GetBlockSize(); site_k++)
                {
                  m++;

                  unsigned int *site_data_p = &map_block_p->site_data[m];

                  // If the site is solid, continue.
                  if (*site_data_p & BIG_NUMBER3)
                  {
                    continue;
                  }

                  int l = -1;

                  switch (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc]))
                  {
                    case FLUID:
                      l = 0;
                      break;
                    case EDGE:
                      l = 1;
                      break;
                    case INLET:
                      l = 2;
                      break;
                    case OUTLET:
                      l = 3;
                      break;
                    case (INLET | EDGE):
                      l = 4;
                      break;
                    case (OUTLET | EDGE):
                      l = 5;
                      break;
                  }

                  ++lSiteIndexOnProc;

                  bool lIsInnerSite = true;

                  // Iterate over all direction vectors.
                  for (unsigned int q = 1; q < D3Q15::NUMVECTORS; q++)
                  {
                    // Find the neighbour site co-ords in this direction.
                    site_t neigh_i = site_i + D3Q15::CX[q];
                    site_t neigh_j = site_j + D3Q15::CY[q];
                    site_t neigh_k = site_k + D3Q15::CZ[q];

                    if (!IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // Find the processor Id for that neighbour.
                    const proc_t
                        * proc_id_p = GetProcIdFromGlobalCoords(util::Vector3D<site_t>(neigh_i,
                                                                                       neigh_j,
                                                                                       neigh_k));

                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                    // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                    // in lbmReadConfig in io.cc.
                    if (proc_id_p == NULL || (int) netTop->GetLocalRank() == *proc_id_p
                        || *proc_id_p == (BIG_NUMBER2))
                    {
                      continue;
                    }

                    lIsInnerSite = false;
                    break;
                  }

                  if (lIsInnerSite)
                  {
                    if (l != 0)
                    {
                      *site_data_p = (unsigned int) (collision_offset[0][l] + innerColsPassed[l]);
                      ++innerColsPassed[l];
                    }
                  }
                  else
                  {
                    *site_data_p = (unsigned int) (collision_offset[1][l] + interColsPassed[l]);
                    ++interColsPassed[l];
                  }
                }
              }
            }
          }
        }
      }

      localLatDat.my_inner_sites = innerSites;

      for (unsigned int m = 0; m < COLLISION_TYPES; ++m)
      {
        localLatDat.my_inter_collisions[m] = interCollisions[m];
        localLatDat.my_inner_collisions[m] = innerCollisions[m];
      }

      localLatDat.SetSharedSiteCount(totalSharedFs);
    }

    void LatticeData::InitialisePointToPointComms(site_t** &lSharedFLocationForEachProc)
    {
      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      // point-to-point communications are performed to match data to be
      // sent to/receive from different partitions; in this way, the
      // communication of the locations of the interface-dependent fluid
      // sites and the identifiers of the distribution functions which
      // propagate to different partitions is avoided (only their values
      // will be communicated). It's here!
      // Allocate the request variable.

      net::Net tempNet;

      for (size_t m = 0; m < neighbouringProcs.size(); m++)
      {
        NeighbouringProcessor *neigh_proc_p = &neighbouringProcs[m];
        // One way send receive.  The lower numbered netTop->ProcessorCount send and the higher numbered ones receive.
        // It seems that, for each pair of processors, the lower numbered one ends up with its own
        // edge sites and directions stored and the higher numbered one ends up with those on the
        // other processor.
        if (neigh_proc_p->Rank > netTop->GetLocalRank())
        {
          tempNet.RequestSend(&lSharedFLocationForEachProc[m][0],
                              neigh_proc_p->SharedFCount * 4,
                              neigh_proc_p->Rank);
        }
        else
        {
          tempNet.RequestReceive(&lSharedFLocationForEachProc[m][0],
                                 neigh_proc_p->SharedFCount * 4,
                                 neigh_proc_p->Rank);
        }
      }

      tempNet.Send();
      tempNet.Receive();
      tempNet.Wait();
    }

    void LatticeData::Initialise()
    {
      // Create a map between the two-level data representation and the 1D
      // compact one is created here.

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      // This rank's site data.
      unsigned int *lThisRankSiteData =
          new unsigned int[globLatDat.fluidSitesOnEachProcessor[netTop->GetLocalRank()]];

      globLatDat.GetThisRankSiteData(lThisRankSiteData);

      // The numbers of inter- and intra-machine neighbouring processors,
      // interface-dependent and independent fluid sites and shared
      // distribution functions of the reference processor are calculated
      // here.  neigh_proc is a static array that is declared in config.h.

      CountCollisionTypes(lThisRankSiteData);

      // the precise interface-dependent data (interface-dependent fluid
      // site locations and identifiers of the distribution functions
      // streamed between different partitions) are collected and the
      // buffers needed for the communications are set from here

      site_t* f_data = new site_t[4 * totalSharedFs];

      // Allocate the index in which to put the distribution functions received from the other
      // process.
      f_recv_iv = new site_t[totalSharedFs];

      site_t** lSharedFLocationForEachProc = new site_t*[neighbouringProcs.size()];

      // Reset to zero again.
      totalSharedFs = 0;

      // Set the remaining neighbouring processor data.
      for (size_t n = 0; n < neighbouringProcs.size(); n++)
      {
        // f_data compacted according to number of shared f_s on each process.
        // f_data will be set later.


        // Site co-ordinates for each of the shared distribution, and the number of
        // the corresponding direction vector. Array is 4 elements for each shared distribution.
        lSharedFLocationForEachProc[n] = &f_data[totalSharedFs << 2];

        // Pointing to a few things, but not setting any variables.
        // FirstSharedF points to start of shared_fs.
        neighbouringProcs[n].FirstSharedF = GetLocalFluidSiteCount() * D3Q15::NUMVECTORS + 1
            + totalSharedFs;

        totalSharedFs += neighbouringProcs[n].SharedFCount;
      }

      neighbourIndexFromProcRank = new proc_t[netTop->GetProcessorCount()];

      for (proc_t m = 0; m < netTop->GetProcessorCount(); m++)
      {
        neighbourIndexFromProcRank[m] = -1;
      }
      // Get neigh_proc_index from ProcessorRankForEachBlockSite.
      for (proc_t m = 0; m < (proc_t) neighbouringProcs.size(); m++)
      {
        neighbourIndexFromProcRank[neighbouringProcs[m].Rank] = m;
      }

      {
        site_t** SharedLocationPerProcByNeighbourId = new site_t*[netTop->GetProcessorCount()];

        for (proc_t ii = 0; ii < netTop->GetProcessorCount(); ++ii)
        {
          if (neighbourIndexFromProcRank[ii] >= 0)
          {
            SharedLocationPerProcByNeighbourId[ii]
                = lSharedFLocationForEachProc[neighbourIndexFromProcRank[ii]];
          }
        }

        InitialiseNeighbourLookup(SharedLocationPerProcByNeighbourId,
                                  netTop->GetLocalRank(),
                                  lThisRankSiteData);

        delete[] SharedLocationPerProcByNeighbourId;
      }

      delete[] lThisRankSiteData;

      InitialisePointToPointComms(lSharedFLocationForEachProc);

      site_t f_count = GetLocalFluidSiteCount() * D3Q15::NUMVECTORS;

      site_t sharedSitesSeen = 0;

      for (size_t m = 0; m < neighbouringProcs.size(); m++)
      {
        NeighbouringProcessor *neigh_proc_p = &neighbouringProcs[m];

        for (site_t n = 0; n < neigh_proc_p->SharedFCount; n++)
        {
          // Get coordinates and direction of the distribution function to be sent to another process.
          site_t* f_data_p = &lSharedFLocationForEachProc[m][n * 4];
          site_t i = f_data_p[0];
          site_t j = f_data_p[1];
          site_t k = f_data_p[2];
          site_t l = f_data_p[3];

          // Correct so that each process has the correct coordinates.
          if (neigh_proc_p->Rank < netTop->GetLocalRank())
          {
            i += D3Q15::CX[l];
            j += D3Q15::CY[l];
            k += D3Q15::CZ[l];
            l = D3Q15::INVERSEDIRECTIONS[l];
          }

          // Get the fluid site number of site that will send data to another process.
          site_t contigSiteId = GetContiguousSiteId(i, j, k);

          // Set f_id to the element in the send buffer that we put the updated
          // distribution functions in.
          localLatDat.SetNeighbourLocation(contigSiteId, (unsigned int) l, ++f_count);

          // Set the place where we put the received distribution functions, which is
          // f_new[number of fluid site that sends, inverse direction].
          f_recv_iv[sharedSitesSeen] = contigSiteId * D3Q15::NUMVECTORS
              + D3Q15::INVERSEDIRECTIONS[l];
          ++sharedSitesSeen;
        }
      }
      // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we
      // are freeing both of those.
      delete[] f_data;

      // Delete the array in which we kept the shared f locations. Don't delete subarrays - these
      // are pointers to elsewhere.
      delete[] lSharedFLocationForEachProc;
    }

    void LatticeData::SendAndReceive(hemelb::net::Net* net)
    {
      for (std::vector<NeighbouringProcessor>::const_iterator it = neighbouringProcs.begin(); it
          != neighbouringProcs.end(); it++)
      {
        // Request the receive into the appropriate bit of FOld.
        net->RequestReceive<distribn_t> (GetFOld( (*it).FirstSharedF),
                                         (int) (*it).SharedFCount,
                                          (*it).Rank);

        // Request the send from the right bit of FNew.
        net->RequestSend<distribn_t> (GetFNew( (*it).FirstSharedF),
                                      (int) (*it).SharedFCount,
                                       (*it).Rank);

      }
    }

    void LatticeData::SwapOldAndNew()
    {
      distribn_t* temp = localLatDat.FOld;
      localLatDat.FOld = localLatDat.FNew;
      localLatDat.FNew = temp;
    }

    void LatticeData::CopyReceived()
    {
      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      for (site_t i = 0; i < totalSharedFs; i++)
      {
        *GetFNew(f_recv_iv[i]) = *GetFOld(neighbouringProcs[0].FirstSharedF + i);
      }
    }

    const site_t* LatticeData::GetFluidSiteCountsOnEachProc() const
    {
      return globLatDat.fluidSitesOnEachProcessor;
    }

    site_t LatticeData::GetFluidSiteCountOnProc(proc_t proc) const
    {
      return globLatDat.fluidSitesOnEachProcessor[proc];
    }

    site_t LatticeData::GetTotalFluidSites() const
    {
      return totalFluidSites;
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteMins() const
    {
      return globLatDat.globalSiteMins;
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteMaxes() const
    {
      return globLatDat.globalSiteMaxes;
    }
  }
}
