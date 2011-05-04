#ifndef HEMELB_GEOMETRY_LATTICEDATA_H
#define HEMELB_GEOMETRY_LATTICEDATA_H

#include <cstdio>

#include "D3Q15.h"
#include "constants.h"
#include "SimConfig.h"
#include "mpiInclude.h"
// TODO Remove the stress type from the data file, so we can remove the dependence
// on LbmParams here.
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace geometry
  {
    class LatticeData
    {
      public:
        enum SiteType
        {
          // These must be consistent with the setup tool
          SOLID_TYPE = 0U,
          FLUID_TYPE = 1U,
          INLET_TYPE = 2U,
          OUTLET_TYPE = 3U
        };

        //TODO Ideally we'd hide implementation details like this.
        // Data about an element of the domain wall
        struct WallData
        {
            // estimated wall normal (if the site is close to the wall);
            double wall_nor[3];
            // cut distances along the 14 non-zero lattice vectors;
            // each one is between 0 and 1 if the surface cuts the corresponding
            // vector or is equal to "NO_VALUE" otherwise
            double cut_dist[D3Q15::NUMVECTORS - 1];
        };

        // Data about each global block in the lattice,
        // site_data[] is an array containing individual lattice site data
        // within a global block.
        struct BlockData
        {
            BlockData()
            {
              ProcessorRankForEachBlockSite = NULL;
              wall_data = NULL;
              site_data = NULL;
            }

            ~BlockData()
            {
              if (ProcessorRankForEachBlockSite != NULL)
              {
                delete[] ProcessorRankForEachBlockSite;
                ProcessorRankForEachBlockSite = NULL;
              }
              if (wall_data != NULL)
              {
                delete[] wall_data;
                wall_data = NULL;
              }
              if (site_data != NULL)
              {
                delete[] site_data;
                site_data = NULL;
              }
            }

            // An array of the ranks on which each lattice site within the block resides.
            proc_t* ProcessorRankForEachBlockSite;
            // Information about wall / inlet / outlet position and orientation for
            // each site.
            WallData *wall_data;
            // The "site data" for each site.
            unsigned int *site_data;
        };

        LatticeData(const bool reserveSteeringCore,
                    site_t* totalFluidSites,
                    site_t siteMins[3],
                    site_t siteMaxes[3],
                    site_t* fluidSitePerProc,
                    lb::LbmParameters* bLbmParams,
                    SimConfig* bSimConfig,
                    double* lReadTime,
                    double* lDecomposeTime);

        void InitialiseNeighbourLookup(site_t** bSharedFLocationForEachProc,
                                       proc_t localRank,
                                       const unsigned int* iSiteDataForThisRank);

        void SetNeighbourLocation(site_t iSiteIndex, unsigned int iDirection, site_t iValue);

        void SetSiteCounts(site_t innerSites,
                           site_t interCollisions[COLLISION_TYPES],
                           site_t innerCollisions[COLLISION_TYPES],
                           site_t sharedSites);

        void SwapOldAndNew();

        const double* GetNormalToWall(site_t iSiteIndex) const;

        site_t GetXSiteCount() const;
        site_t GetYSiteCount() const;
        site_t GetZSiteCount() const;

        site_t GetXBlockCount() const;
        site_t GetYBlockCount() const;
        site_t GetZBlockCount() const;

        unsigned int GetLog2BlockSize() const;

        site_t GetBlockSize() const;
        site_t GetBlockCount() const;

        site_t GetSitesPerBlockVolumeUnit() const;

        site_t GetBlockIdFromBlockCoords(site_t i, site_t j, site_t k) const;

        bool IsValidLatticeSite(site_t i, site_t j, site_t k) const;

        const proc_t* GetProcIdFromGlobalCoords(site_t siteI, site_t siteJ, site_t siteK) const;

        BlockData* GetBlock(site_t blockNumber) const;

        distribn_t* GetFOld(site_t siteNumber) const;
        distribn_t* GetFNew(site_t siteNumber) const;
        site_t GetLocalFluidSiteCount() const;
        SiteType GetSiteType(site_t iSiteIndex) const;
        int GetBoundaryId(site_t iSiteIndex) const;
        site_t GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const;
        bool HasBoundary(site_t iSiteIndex, int iDirection) const;
        double GetCutDistance(site_t iSiteIndex, int iDirection) const;
        unsigned int GetSiteData(site_t iSiteIndex) const;
        unsigned int GetContiguousSiteId(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;
        site_t GetInnerSiteCount() const;
        site_t GetInnerCollisionCount(unsigned int collisionType) const;
        site_t GetInterCollisionCount(unsigned int collisionType) const;
        unsigned int GetCollisionType(unsigned int site_data) const;

      private:
        class LocalLatticeData
        {
          public:
            LocalLatticeData();
            ~LocalLatticeData();

            void Initialise(site_t iLocalFluidSites);

            site_t GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const;
            double GetCutDistance(site_t iSiteIndex, int iDirection) const;
            bool HasBoundary(site_t iSiteIndex, int iDirection) const;
            int GetBoundaryId(site_t iSiteIndex) const;
            const double *GetNormalToWall(site_t iSiteIndex) const;
            SiteType GetSiteType(site_t iSiteIndex) const;
            site_t GetLocalFluidSiteCount() const;

            void SetNeighbourLocation(site_t iSiteIndex, unsigned int iDirection, site_t iValue);
            void SetWallNormal(site_t iSiteIndex, const double iNormal[3]);
            void
            SetDistanceToWall(site_t iSiteIndex, const double iCutDistance[D3Q15::NUMVECTORS - 1]);

            void SetSharedSiteCount(site_t iSharedCount);

          public:
            site_t my_inner_sites;
            site_t my_inner_collisions[COLLISION_TYPES];
            site_t my_inter_collisions[COLLISION_TYPES];

            distribn_t *FOld;
            distribn_t *FNew;

            // TODO sadly this has to be public, due to some budgetry in the way we determine site type.
            // SiteType || FluidSite and SiteType && FluidSite have different significances...
            unsigned int *mSiteData;

          private:
            site_t LocalFluidSites;
            site_t* mFNeighbours;
            double *mDistanceToWall;
            double *mWallNormalAtSite;
        };

        class GlobalLatticeData
        {
            friend class BlockCounter;

          public:
            void SetBasicDetails(site_t iBlocksX,
                                 site_t iBlocksY,
                                 site_t iBlocksZ,
                                 site_t iBlockSize);

            site_t GetXSiteCount() const;
            site_t GetYSiteCount() const;
            site_t GetZSiteCount() const;
            site_t GetXBlockCount() const;
            site_t GetYBlockCount() const;
            site_t GetZBlockCount() const;

            site_t GetBlockSize() const;
            site_t GetBlockCount() const;

            site_t GetSitesPerBlockVolumeUnit() const;

            bool IsValidLatticeSite(site_t i, site_t j, site_t k) const;

            BlockData * Blocks;

            ~GlobalLatticeData();

            // Returns the type of collision/streaming update for the fluid site
            // with data "site_data".
            unsigned int GetCollisionType(unsigned int site_data) const;

            // Function that finds the pointer to the rank on which a particular site
            // resides. If the site is in an empty block, return NULL.
            const proc_t
            * GetProcIdFromGlobalCoords(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;

            // Function that gets the index of a block from its coordinates.
            site_t GetBlockIdFromBlockCoords(site_t blockI, site_t blockJ, site_t blockK) const;

            unsigned int GetSiteData(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;

          public:
            // TODO public temporarily, until all usages are internal to the class.
            unsigned int Log2BlockSize;

          private:
            site_t mSitesPerBlockVolumeUnit;
            site_t mBlockCount;
            site_t mSitesX, mSitesY, mSitesZ;
            site_t mBlocksX, mBlocksY, mBlocksZ;
            site_t mBlockSize;
        };

        class BlockCounter
        {
          public:
            BlockCounter(const GlobalLatticeData* iGlobLatDat, site_t iStartNumber)
            {
              mBlockNumber = iStartNumber;
              mGlobLatDat = iGlobLatDat;
            }

            void operator++()
            {
              mBlockNumber++;
            }

            void operator++(int in)
            {
              mBlockNumber++;
            }

            operator site_t()
            {
              return mBlockNumber;
            }

            bool operator<(site_t iUpperLimit) const
            {
              return mBlockNumber < iUpperLimit;
            }

            site_t GetICoord()
            {
              return (mBlockNumber - (mBlockNumber % (mGlobLatDat->GetYBlockCount()
                  * mGlobLatDat->GetZBlockCount()))) / (mGlobLatDat->GetYBlockCount()
                  * mGlobLatDat->GetZBlockCount());
            }

            site_t GetJCoord()
            {
              site_t lTemp = mBlockNumber % (mGlobLatDat->GetYBlockCount()
                  * mGlobLatDat->GetZBlockCount());
              return (lTemp - (lTemp % mGlobLatDat->GetZBlockCount()))
                  / mGlobLatDat->GetZBlockCount();
            }

            site_t GetKCoord()
            {
              return mBlockNumber % mGlobLatDat->GetZBlockCount();
            }

            site_t GetICoord(site_t iSiteI)
            {
              return (GetICoord() << mGlobLatDat->Log2BlockSize) + iSiteI;
            }

            site_t GetJCoord(site_t iSiteJ)
            {
              return (GetJCoord() << mGlobLatDat->Log2BlockSize) + iSiteJ;
            }

            site_t GetKCoord(site_t iSiteK)
            {
              return (GetKCoord() << mGlobLatDat->Log2BlockSize) + iSiteK;
            }

          private:
            const GlobalLatticeData* mGlobLatDat;
            site_t mBlockNumber;
        };

        class GeometryReader
        {
          public:
            GeometryReader(const bool reserveSteeringCore);
            ~GeometryReader();

            void LoadAndDecompose(GlobalLatticeData* bGlobalLatticeData,
                                  site_t* totalFluidSites,
                                  site_t siteMins[3],
                                  site_t siteMaxes[3],
                                  site_t* fluidSitePerProc,
                                  lb::LbmParameters* bLbmParams,
                                  SimConfig* bSimConfig,
                                  double* lReadTime,
                                  double* lDecomposeTime);

          private:
            struct BlockLocation
            {
                site_t i, j, k;
            };

            void ReadPreamble(MPI_File xiFile,
                              lb::LbmParameters* bParams,
                              GlobalLatticeData* bGlobalLatticeData);

            void ReadHeader(MPI_File xiFile,
                            site_t iBlockCount,
                            site_t* sitesInEachBlock,
                            unsigned int* bytesUsedByBlockInDataFile);

            void BlockDecomposition(const site_t iBlockCount,
                                    const GlobalLatticeData* iGlobLatDat,
                                    const site_t* fluidSitePerBlock,
                                    proc_t* initialProcForEachBlock);

            void DivideBlocks(site_t unassignedBlocks,
                              site_t totalBlockCount,
                              proc_t unitCount,
                              site_t* blocksOnEachUnit,
                              proc_t* unitForEachBlock,
                              const site_t* fluidSitesPerBlock,
                              const GlobalLatticeData* iGlobLatDat);

            void ReadInLocalBlocks(MPI_File iFile,
                                   const unsigned int* bytesPerBlock,
                                   const proc_t* unitForEachBlock,
                                   const proc_t localRank,
                                   const GlobalLatticeData* iGlobLatDat);

            void OptimiseDomainDecomposition(const site_t* sitesPerBlock,
                                             const unsigned int* bytesPerBlock,
                                             const proc_t* procForEachBlock,
                                             MPI_File iFile,
                                             GlobalLatticeData* bGlobLatDat);

            MPI_Comm mTopologyComm;
            MPI_Group mTopologyGroup;
            int mTopologyRank;
            unsigned int mTopologySize;
            int mGlobalRank;
            bool mParticipateInTopology;
        };

        void SetSiteData(site_t siteIndex, unsigned int siteData);
        void SetWallNormal(site_t siteIndex, double normal[3]);
        void SetWallDistance(site_t siteIndex, double cutDistance[D3Q15::NUMVECTORS - 1]);

        LocalLatticeData localLatDat;
        GlobalLatticeData globLatDat;
    };
  }
}

#endif /* HEMELB_GEOMETRY_LATTICEDATA_H */
