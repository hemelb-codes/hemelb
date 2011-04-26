#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace geometry
  {
    void LatticeData::GlobalLatticeData::SetBasicDetails(unsigned int iBlocksX,
                                                         unsigned int iBlocksY,
                                                         unsigned int iBlocksZ,
                                                         unsigned int iBlockSize)
    {
      mBlocksX = iBlocksX;
      mBlocksY = iBlocksY;
      mBlocksZ = iBlocksZ;
      mBlockSize = iBlockSize;

      mSitesX = mBlocksX * mBlockSize;
      mSitesY = mBlocksY * mBlockSize;
      mSitesZ = mBlocksZ * mBlockSize;

      mSitesPerBlockVolumeUnit = mBlockSize * mBlockSize * mBlockSize;

      // A shift value we'll need later = log_2(block_size)
      unsigned int i = mBlockSize;
      Log2BlockSize = 0;
      while (i > 1)
      {
        i >>= 1;
        ++Log2BlockSize;
      }

      mBlockCount = mBlocksX * mBlocksY * mBlocksZ;

      Blocks = new BlockData[mBlockCount];
    }

    unsigned int LatticeData::GlobalLatticeData::GetXSiteCount() const
    {
      return mSitesX;
    }

    unsigned int LatticeData::GlobalLatticeData::GetYSiteCount() const
    {
      return mSitesY;
    }

    unsigned int LatticeData::GlobalLatticeData::GetZSiteCount() const
    {
      return mSitesZ;
    }

    unsigned int LatticeData::GlobalLatticeData::GetXBlockCount() const
    {
      return mBlocksX;
    }

    unsigned int LatticeData::GlobalLatticeData::GetYBlockCount() const
    {
      return mBlocksY;
    }

    unsigned int LatticeData::GlobalLatticeData::GetZBlockCount() const
    {
      return mBlocksZ;
    }

    unsigned int LatticeData::GlobalLatticeData::GetBlockSize() const
    {
      return mBlockSize;
    }

    unsigned int LatticeData::GlobalLatticeData::GetBlockCount() const
    {
      return mBlockCount;
    }

    unsigned int LatticeData::GlobalLatticeData::GetSitesPerBlockVolumeUnit() const
    {
      return mSitesPerBlockVolumeUnit;
    }

    bool LatticeData::GlobalLatticeData::IsValidLatticeSite(unsigned int i,
                                                            unsigned int j,
                                                            unsigned int k) const
    {
      return i < mSitesX && j < mSitesY && k < mSitesZ;
    }

    LatticeData::GlobalLatticeData::~GlobalLatticeData()
    {
      delete[] Blocks;
    }

    // Returns the type of collision/streaming update for the fluid site
    // with data "site_data".
    unsigned int LatticeData::GlobalLatticeData::GetCollisionType(unsigned int site_data) const
    {
      unsigned int boundary_type;

      if (site_data == LatticeData::FLUID_TYPE)
      {
        return FLUID;
      }
      boundary_type = site_data & SITE_TYPE_MASK;

      if (boundary_type == LatticeData::FLUID_TYPE)
      {
        return EDGE;
      }
      if (! (site_data & PRESSURE_EDGE_MASK))
      {
        if (boundary_type == LatticeData::INLET_TYPE)
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
        if (boundary_type == LatticeData::INLET_TYPE)
        {
          return INLET | EDGE;
        }
        else
        {
          return OUTLET | EDGE;
        }
      }
    }

    // Function that finds the rank on which a particular site
    // resides. If the site is in an empty block, return -1.
    int LatticeData::GlobalLatticeData::GetProcIdFromGlobalCoords(unsigned int iSiteI,
                                                                  unsigned int iSiteJ,
                                                                  unsigned int iSiteK) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      unsigned int i = iSiteI >> Log2BlockSize;
      unsigned int j = iSiteJ >> Log2BlockSize;
      unsigned int k = iSiteK >> Log2BlockSize;

      // Get the block from the block identifiers.
      BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(i, j, k)];

      // If an empty (solid) block is addressed, return a NULL pointer.
      if (lBlock->ProcessorRankForEachBlockSite == NULL)
      {
        return -1;
      }
      else
      {
        // Find site coordinates within the block
        unsigned int ii = iSiteI - (i << Log2BlockSize);
        unsigned int jj = iSiteJ - (j << Log2BlockSize);
        unsigned int kk = iSiteK - (k << Log2BlockSize);

        // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
        // mProcessorsForEachBlock)
        return lBlock->ProcessorRankForEachBlockSite[ ( ( (ii << Log2BlockSize) + jj)
            << Log2BlockSize) + kk];
      }
    }

    // Function that gets the index of a block from its coordinates.
    unsigned int LatticeData::GlobalLatticeData::GetBlockIdFromBlockCoords(unsigned int blockI,
                                                                           unsigned int blockJ,
                                                                           unsigned int blockK) const
    {
      // Get the block from the block identifiers.
      return (blockI * mBlocksY + blockJ) * mBlocksZ + blockK;
    }

    // Function to get a pointer to the site_data for a site.
    // If the site is in an empty block, return NULL.

    unsigned int LatticeData::GlobalLatticeData::GetSiteData(unsigned int iSiteI,
                                                             unsigned int iSiteJ,
                                                             unsigned int iSiteK) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      unsigned int i = iSiteI >> Log2BlockSize;
      unsigned int j = iSiteJ >> Log2BlockSize;
      unsigned int k = iSiteK >> Log2BlockSize;

      // Pointer to the block
      BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(i, j, k)];

      // Find site coordinates within the block
      unsigned int ii = iSiteI - (i << Log2BlockSize);
      unsigned int jj = iSiteJ - (j << Log2BlockSize);
      unsigned int kk = iSiteK - (k << Log2BlockSize);

      // Return pointer to site_data[site]
      return lBlock->site_data[ ( ( (ii << Log2BlockSize) + jj) << Log2BlockSize) + kk];
    }
  }
}
