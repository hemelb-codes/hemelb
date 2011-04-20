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
            // vector or is equal to "BIG_NUMBER" otherwise
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
            int *ProcessorRankForEachBlockSite;
            // Information about wall / inlet / outlet position and orientation for
            // each site.
            WallData *wall_data;
            // The "site data" for each site.
            unsigned int *site_data;
        };

        LatticeData(const bool reserveSteeringCore,
                    int* totalFluidSites,
                    unsigned int siteMins[3],
                    unsigned int siteMaxes[3],
                    unsigned int* fluidSitePerProc,
                    lb::LbmParameters* bLbmParams,
                    SimConfig* bSimConfig,
                    double* lReadTime,
                    double* lDecomposeTime);

        void InitialiseNeighbourLookup(int ** bSharedFLocationForEachProc,
                                       int localRank,
                                       const unsigned int* iSiteDataForThisRank);

        void SetNeighbourLocation(unsigned int iSiteIndex,
                                  unsigned int iDirection,
                                  unsigned int iValue);

        void SetSiteCounts(unsigned int innerSites,
                           unsigned int interCollisions[COLLISION_TYPES],
                           unsigned int innerCollisions[COLLISION_TYPES],
                           unsigned int sharedSites);

        void SwapOldAndNew();

        const double* GetNormalToWall(int iSiteIndex) const;

        unsigned int GetXSiteCount() const;
        unsigned int GetYSiteCount() const;
        unsigned int GetZSiteCount() const;

        unsigned int GetXBlockCount() const;
        unsigned int GetYBlockCount() const;
        unsigned int GetZBlockCount() const;

        unsigned int GetLog2BlockSize() const;

        unsigned int GetBlockSize() const;
        unsigned int GetBlockCount() const;

        unsigned int GetSitesPerBlockVolumeUnit() const;

        unsigned int GetBlockIdFromBlockCoords(unsigned int i, unsigned int j, unsigned int k) const;

        bool IsValidLatticeSite(unsigned int i, unsigned int j, unsigned int k) const;

        int*
        GetProcIdFromGlobalCoords(unsigned int siteI, unsigned int siteJ, unsigned int siteK) const;

        BlockData* GetBlock(unsigned int blockNumber) const;

        double* GetFOld(unsigned int siteNumber) const;
        double* GetFNew(unsigned int siteNumber) const;
        unsigned int GetLocalFluidSiteCount() const;
        SiteType GetSiteType(int iSiteIndex) const;
        int GetBoundaryId(int iSiteIndex) const;
        int GetStreamedIndex(int iSiteIndex, int iDirectionIndex) const;
        bool HasBoundary(int iSiteIndex, int iDirection) const;
        double GetCutDistance(int iSiteIndex, int iDirection) const;
        unsigned int GetSiteData(unsigned int iSiteIndex) const;
        unsigned int GetContiguousSiteId(unsigned int iSiteI,
                                         unsigned int iSiteJ,
                                         unsigned int iSiteK) const;
        unsigned int GetInnerSiteCount() const;
        unsigned int GetInnerCollisionCount(unsigned int collisionType) const;
        unsigned int GetInterCollisionCount(unsigned int collisionType) const;
        unsigned int GetCollisionType(unsigned int site_data) const;

      private:
        class LocalLatticeData
        {
          public:
            LocalLatticeData();
            ~LocalLatticeData();

            void Initialise(unsigned int iLocalFluidSites);

            int GetStreamedIndex(int iSiteIndex, int iDirectionIndex) const;
            double GetCutDistance(int iSiteIndex, int iDirection) const;
            bool HasBoundary(int iSiteIndex, int iDirection) const;
            int GetBoundaryId(int iSiteIndex) const;
            const double *GetNormalToWall(int iSiteIndex) const;
            SiteType GetSiteType(int iSiteIndex) const;
            unsigned int GetLocalFluidSiteCount() const;

            void SetNeighbourLocation(unsigned int iSiteIndex,
                                      unsigned int iDirection,
                                      unsigned int iValue);
            void SetWallNormal(int iSiteIndex, const double iNormal[3]);
            void
            SetDistanceToWall(int iSiteIndex, const double iCutDistance[D3Q15::NUMVECTORS - 1]);

            void SetSharedSiteCount(int iSharedCount);

          public:
            unsigned int my_inner_sites;
            unsigned int my_inner_collisions[COLLISION_TYPES];
            unsigned int my_inter_collisions[COLLISION_TYPES];

            double *FOld;
            double *FNew;

            // TODO sadly this has to be public, due to some budgetry in the way we determine site type.
            // SiteType || FluidSite and SiteType && FluidSite have different significances...
            unsigned int *mSiteData;

          private:
            unsigned int LocalFluidSites;
            unsigned int *mFNeighbours;
            double *mDistanceToWall;
            double *mWallNormalAtSite;
        };

        class GlobalLatticeData
        {
            friend class BlockCounter;

          public:
            void SetBasicDetails(unsigned int iBlocksX,
                                 unsigned int iBlocksY,
                                 unsigned int iBlocksZ,
                                 unsigned int iBlockSize);

            unsigned int GetXSiteCount() const;
            unsigned int GetYSiteCount() const;
            unsigned int GetZSiteCount() const;
            unsigned int GetXBlockCount() const;
            unsigned int GetYBlockCount() const;
            unsigned int GetZBlockCount() const;

            unsigned int GetBlockSize() const;
            unsigned int GetBlockCount() const;

            unsigned int GetSitesPerBlockVolumeUnit() const;

            bool IsValidLatticeSite(unsigned int i, unsigned int j, unsigned int k) const;

            BlockData * Blocks;

            ~GlobalLatticeData();

            // Returns the type of collision/streaming update for the fluid site
            // with data "site_data".
            unsigned int GetCollisionType(unsigned int site_data) const;

            // Function that finds the pointer to the rank on which a particular site
            // resides. If the site is in an empty block, return NULL.
            int * GetProcIdFromGlobalCoords(unsigned int iSiteI,
                                            unsigned int iSiteJ,
                                            unsigned int iSiteK) const;

            // Function that gets the index of a block from its coordinates.
            unsigned int GetBlockIdFromBlockCoords(unsigned int blockI,
                                                   unsigned int blockJ,
                                                   unsigned int blockK) const;

            unsigned int
            GetSiteData(unsigned int iSiteI, unsigned int iSiteJ, unsigned int iSiteK) const;
          public:
            // TODO public temporarily, until all usages are internal to the class.
            unsigned int Log2BlockSize;

          private:
            unsigned int mSitesPerBlockVolumeUnit;
            unsigned int mBlockCount;
            unsigned int mSitesX, mSitesY, mSitesZ;
            unsigned int mBlocksX, mBlocksY, mBlocksZ;
            unsigned int mBlockSize;
        };

        class BlockCounter
        {
          public:
            BlockCounter(const GlobalLatticeData* iGlobLatDat, unsigned int iStartNumber)
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

            operator int()
            {
              return mBlockNumber;
            }

            bool operator<(unsigned int iUpperLimit) const
            {
              return mBlockNumber < iUpperLimit;
            }

            int GetICoord()
            {
              return (mBlockNumber - (mBlockNumber % (mGlobLatDat->GetYBlockCount()
                  * mGlobLatDat->GetZBlockCount()))) / (mGlobLatDat->GetYBlockCount()
                  * mGlobLatDat->GetZBlockCount());
            }

            int GetJCoord()
            {
              int lTemp = mBlockNumber % (mGlobLatDat->GetYBlockCount()
                  * mGlobLatDat->GetZBlockCount());
              return (lTemp - (lTemp % mGlobLatDat->GetZBlockCount()))
                  / mGlobLatDat->GetZBlockCount();
            }

            int GetKCoord()
            {
              return mBlockNumber % mGlobLatDat->GetZBlockCount();
            }

            int GetICoord(int iSiteI)
            {
              return (GetICoord() << mGlobLatDat->Log2BlockSize) + iSiteI;
            }

            int GetJCoord(int iSiteJ)
            {
              return (GetJCoord() << mGlobLatDat->Log2BlockSize) + iSiteJ;
            }

            int GetKCoord(int iSiteK)
            {
              return (GetKCoord() << mGlobLatDat->Log2BlockSize) + iSiteK;
            }

          private:
            const GlobalLatticeData* mGlobLatDat;
            unsigned int mBlockNumber;
        };

        class GeometryReader
        {
          public:
            GeometryReader(const bool reserveSteeringCore);
            ~GeometryReader();

            void LoadAndDecompose(GlobalLatticeData* bGlobalLatticeData,
                                  int* totalFluidSites,
                                  unsigned int siteMins[3],
                                  unsigned int siteMaxes[3],
                                  unsigned int* fluidSitePerProc,
                                  lb::LbmParameters* bLbmParams,
                                  SimConfig* bSimConfig,
                                  double* lReadTime,
                                  double* lDecomposeTime);

          private:
            struct BlockLocation
            {
                int i, j, k;
            };

            void ReadPreamble(MPI_File xiFile,
                              lb::LbmParameters* bParams,
                              GlobalLatticeData* bGlobalLatticeData);

            void ReadHeader(MPI_File xiFile,
                            unsigned int iBlockCount,
                            unsigned int* sitesInEachBlock,
                            unsigned int* bytesUsedByBlockInDataFile);

            void BlockDecomposition(const unsigned int iBlockCount,
                                    const GlobalLatticeData* iGlobLatDat,
                                    const unsigned int* fluidSitePerBlock,
                                    int* initialProcForEachBlock);

            void DivideBlocks(unsigned int unassignedBlocks,
                              unsigned int totalBlockCount,
                              unsigned int unitCount,
                              unsigned int* blocksOnEachUnit,
                              int* unitForEachBlock,
                              const unsigned int* fluidSitesPerBlock,
                              const GlobalLatticeData* iGlobLatDat);

            void ReadInLocalBlocks(MPI_File iFile,
                                   const unsigned int* bytesPerBlock,
                                   const int* unitForEachBlock,
                                   const unsigned int localRank,
                                   const GlobalLatticeData* iGlobLatDat);

            void OptimiseDomainDecomposition(const unsigned int* sitesPerBlock,
                                             const unsigned int* bytesPerBlock,
                                             const int* procForEachBlock,
                                             MPI_File iFile,
                                             lb::LbmParameters* bLbmParams,
                                             GlobalLatticeData* bGlobLatDat);

            MPI_Comm mTopologyComm;
            MPI_Group mTopologyGroup;
            int mTopologyRank;
            unsigned int mTopologySize;
            int mGlobalRank;
        };

        void SetSiteData(unsigned int siteIndex, unsigned int siteData);
        void SetWallNormal(unsigned int siteIndex, double normal[3]);
        void SetWallDistance(unsigned int siteIndex, double cutDistance[D3Q15::NUMVECTORS - 1]);

        LocalLatticeData localLatDat;
        GlobalLatticeData globLatDat;
    };
  }
}

#endif /* HEMELB_GEOMETRY_LATTICEDATA_H */
