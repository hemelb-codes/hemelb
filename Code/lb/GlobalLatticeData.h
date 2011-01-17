#ifndef HEMELB_LB_GLOBALLATTICEDATA_H
#define HEMELB_LB_GLOBALLATTICEDATA_H

#include "D3Q15.h"
#include "lb/LocalLatticeData.h"
#include <cstdlib>

namespace hemelb
{
  namespace lb
  {
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
            delete[] ProcessorRankForEachBlockSite;
          if (wall_data != NULL)
            delete[] wall_data;
          if (site_data != NULL)
            delete[] site_data;
        }

        // An array of the ranks on which each lattice site within the block resides.
        int *ProcessorRankForEachBlockSite;
        // Information about wall / inlet / outlet position and orientation for
        // each site.
        WallData *wall_data;
        // The "site data" for each site.
        unsigned int *site_data;
    };

    class GlobalLatticeData
    {
      public:
        void SetBasicDetails(int iBlocksX,
                             int iBlocksY,
                             int iBlocksZ,
                             int iBlockSize)
        {
          mBlocksX = iBlocksX;
          mBlocksY = iBlocksY;
          mBlocksZ = iBlocksZ;
          mBlockSize = iBlockSize;

          mSitesX = mBlocksX * mBlockSize;
          mSitesY = mBlocksY * mBlockSize;
          mSitesZ = mBlocksZ * mBlockSize;

          SitesPerBlockVolumeUnit = mBlockSize * mBlockSize * mBlockSize;

          // A shift value we'll need later = log_2(block_size)
          int i = mBlockSize;
          Log2BlockSize = 0;
          while (i > 1)
          {
            i >>= 1;
            ++Log2BlockSize;
          }

          mBlockCount = mBlocksX * mBlocksY * mBlocksZ;

          Blocks = new BlockData[mBlockCount];
        }

        int GetXSiteCount() const
        {
          return mSitesX;
        }

        int GetYSiteCount() const
        {
          return mSitesY;
        }

        int GetZSiteCount() const
        {
          return mSitesZ;
        }

        int GetXBlockCount() const
        {
          return mBlocksX;
        }

        int GetYBlockCount() const
        {
          return mBlocksY;
        }

        int GetZBlockCount() const
        {
          return mBlocksZ;
        }

        int GetBlockSize() const
        {
          return mBlockSize;
        }

        int GetBlockCount() const
        {
          return mBlockCount;
        }

        BlockData * Blocks;

        ~GlobalLatticeData()
        {
          delete[] Blocks;
        }

        // Returns the type of collision/streaming update for the fluid site
        // with data "site_data".
        unsigned int GetCollisionType(unsigned int site_data) const
        {
          unsigned int boundary_type;

          if (site_data == hemelb::lb::FLUID_TYPE)
          {
            return FLUID;
          }
          boundary_type = site_data & SITE_TYPE_MASK;

          if (boundary_type == hemelb::lb::FLUID_TYPE)
          {
            return EDGE;
          }
          if (! (site_data & PRESSURE_EDGE_MASK))
          {
            if (boundary_type == hemelb::lb::INLET_TYPE)
            {
              return INLET;
            }
            else
            {
              return OUTLET;
            }
          }
          else
          {
            if (boundary_type == hemelb::lb::INLET_TYPE)
            {
              return INLET | EDGE;
            }
            else
            {
              return OUTLET | EDGE;
            }
          }
        }

        // Function that finds the pointer to the rank on which a particular site
        // resides. If the site is in an empty block, return NULL.
        int * GetProcIdFromGlobalCoords(int iSiteI, int iSiteJ, int iSiteK) const
        {
          // If the given site location is outside the bounding box return a NULL
          // pointer.
          if (iSiteI < 0 || iSiteI >= mSitesX || iSiteJ < 0 || iSiteJ
              >= mSitesY || iSiteK < 0 || iSiteK >= mSitesZ)
          {
            return NULL;
          }

          // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
          int i = iSiteI >> Log2BlockSize;
          int j = iSiteJ >> Log2BlockSize;
          int k = iSiteK >> Log2BlockSize;

          // Get the block from the block identifiers.
          BlockData * lBlock = &Blocks[ (i * mBlocksY + j) * mBlocksZ + k];

          // If an empty (solid) block is addressed, return a NULL pointer.
          if (lBlock->ProcessorRankForEachBlockSite == NULL)
          {
            return NULL;
          }
          else
          {
            // Find site coordinates within the block
            int ii = iSiteI - (i << Log2BlockSize);
            int jj = iSiteJ - (j << Log2BlockSize);
            int kk = iSiteK - (k << Log2BlockSize);

            // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
            // mProcessorsForEachBlock)
            return &lBlock->ProcessorRankForEachBlockSite[ ( ( (ii
                << Log2BlockSize) + jj) << Log2BlockSize) + kk];
          }
        }

        // Function to get a pointer to the site_data for a site.
        // If the site is in an empty block, return NULL.

        const unsigned int * GetSiteData(int iSiteI, int iSiteJ, int iSiteK) const
        {
          // If site is out of the bounding box, return NULL.
          if (iSiteI < 0 || iSiteI >= mSitesX || iSiteJ < 0 || iSiteJ
              >= mSitesY || iSiteK < 0 || iSiteK >= mSitesZ)
          {
            return NULL;
          }

          // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
          int i = iSiteI >> Log2BlockSize;
          int j = iSiteJ >> Log2BlockSize;
          int k = iSiteK >> Log2BlockSize;

          // Pointer to the block
          BlockData * lBlock = &Blocks[ (i * mBlocksY + j) * mBlocksZ + k];

          // if an empty (solid) block is addressed
          if (lBlock->site_data == NULL)
          {
            return NULL;
          }
          else
          {
            // Find site coordinates within the block
            int ii = iSiteI - (i << Log2BlockSize);
            int jj = iSiteJ - (j << Log2BlockSize);
            int kk = iSiteK - (k << Log2BlockSize);

            // Return pointer to site_data[site]
            return &lBlock->site_data[ ( ( (ii << Log2BlockSize) + jj)
                << Log2BlockSize) + kk];
          }
        }

      public:
        // TODO public temporarily, until all usages are internal to the class.
        int Log2BlockSize;
        int SitesPerBlockVolumeUnit;

      private:
        int mBlockCount;
        int mSitesX, mSitesY, mSitesZ;
        int mBlocksX, mBlocksY, mBlocksZ;
        int mBlockSize;
    };
  }
}

#endif /* HEMELB_LB_GLOBALLATTICEDATA_H */
