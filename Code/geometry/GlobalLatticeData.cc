#include "geometry/LatticeData.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace geometry
  {
    void LatticeData::GlobalLatticeData::SetBasicDetails(site_t iBlocksX,
                                                         site_t iBlocksY,
                                                         site_t iBlocksZ,
                                                         site_t iBlockSize)
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
      site_t i = mBlockSize;
      Log2BlockSize = 0;
      while (i > 1)
      {
        i >>= 1;
        ++Log2BlockSize;
      }

      mBlockCount = mBlocksX * mBlocksY * mBlocksZ;

      Blocks = new BlockData[mBlockCount];
    }

    site_t LatticeData::GlobalLatticeData::GetXSiteCount() const
    {
      return mSitesX;
    }

    site_t LatticeData::GlobalLatticeData::GetYSiteCount() const
    {
      return mSitesY;
    }

    site_t LatticeData::GlobalLatticeData::GetZSiteCount() const
    {
      return mSitesZ;
    }

    site_t LatticeData::GlobalLatticeData::GetXBlockCount() const
    {
      return mBlocksX;
    }

    site_t LatticeData::GlobalLatticeData::GetYBlockCount() const
    {
      return mBlocksY;
    }

    site_t LatticeData::GlobalLatticeData::GetZBlockCount() const
    {
      return mBlocksZ;
    }

    site_t LatticeData::GlobalLatticeData::GetBlockSize() const
    {
      return mBlockSize;
    }

    site_t LatticeData::GlobalLatticeData::GetBlockCount() const
    {
      return mBlockCount;
    }

    site_t LatticeData::GlobalLatticeData::GetSitesPerBlockVolumeUnit() const
    {
      return mSitesPerBlockVolumeUnit;
    }

    bool LatticeData::GlobalLatticeData::IsValidLatticeSite(site_t i, site_t j, site_t k) const
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

    // Function that finds the pointer to the rank on which a particular site
    // resides. If the site is in an empty block, return NULL.
    const proc_t* LatticeData::GlobalLatticeData::GetProcIdFromGlobalCoords(site_t iSiteI,
                                                                            site_t iSiteJ,
                                                                            site_t iSiteK) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      site_t i = iSiteI >> Log2BlockSize;
      site_t j = iSiteJ >> Log2BlockSize;
      site_t k = iSiteK >> Log2BlockSize;

      // Get the block from the block identifiers.
      BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(i, j, k)];

      // If an empty (solid) block is addressed, return a NULL pointer.
      if (lBlock->ProcessorRankForEachBlockSite == NULL)
      {
        return NULL;
      }
      else
      {
        // Find site coordinates within the block
        site_t ii = iSiteI - (i << Log2BlockSize);
        site_t jj = iSiteJ - (j << Log2BlockSize);
        site_t kk = iSiteK - (k << Log2BlockSize);

        // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
        // mProcessorsForEachBlock)
        return &lBlock->ProcessorRankForEachBlockSite[ ( ( (ii << Log2BlockSize) + jj)
            << Log2BlockSize) + kk];
      }
    }

    // Function that gets the index of a block from its coordinates.
    site_t LatticeData::GlobalLatticeData::GetBlockIdFromBlockCoords(site_t blockI,
                                                                     site_t blockJ,
                                                                     site_t blockK) const
    {
      // Get the block from the block identifiers.
      return (blockI * mBlocksY + blockJ) * mBlocksZ + blockK;
    }

    // Function to get a pointer to the site_data for a site.
    // If the site is in an empty block, return NULL.

    void LatticeData::GlobalLatticeData::GetBlockIJK(site_t block, site_t* i, site_t* j, site_t* k) const
    {
      *k = block % GetZBlockCount();
      site_t ij = block / GetZBlockCount();
      *j = ij % GetYBlockCount();
      *i = ij / GetYBlockCount();
    }

    site_t LatticeData::GlobalLatticeData::GetSiteCoord(site_t block, site_t localSiteCoord) const
    {
      return (block << Log2BlockSize) + localSiteCoord;
    }

    unsigned int LatticeData::GlobalLatticeData::GetSiteData(site_t iSiteI,
                                                             site_t iSiteJ,
                                                             site_t iSiteK) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      site_t i = iSiteI >> Log2BlockSize;
      site_t j = iSiteJ >> Log2BlockSize;
      site_t k = iSiteK >> Log2BlockSize;

      // Pointer to the block
      BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(i, j, k)];

      // Find site coordinates within the block
      site_t ii = iSiteI - (i << Log2BlockSize);
      site_t jj = iSiteJ - (j << Log2BlockSize);
      site_t kk = iSiteK - (k << Log2BlockSize);

      // Return pointer to site_data[site]
      return lBlock->site_data[ ( ( (ii << Log2BlockSize) + jj) << Log2BlockSize) + kk];
    }

    void LatticeData::GlobalLatticeData::ReadBlock(const site_t block, io::XdrReader* reader)
    {
      if (Blocks[block].site_data == NULL)
      {
        Blocks[block].site_data = new unsigned int[GetSitesPerBlockVolumeUnit()];
      }
      if (Blocks[block].ProcessorRankForEachBlockSite == NULL)
      {
        Blocks[block].ProcessorRankForEachBlockSite = new proc_t[GetSitesPerBlockVolumeUnit()];
      }

      site_t m = -1;
      site_t blockI, blockJ, blockK;

      GetBlockIJK(block, &blockI, &blockJ, &blockK);

      for (site_t ii = 0; ii < GetBlockSize(); ii++)
      {
        site_t site_i = GetSiteCoord(block, ii);

        for (site_t jj = 0; jj < GetBlockSize(); jj++)
        {
          site_t site_j = GetSiteCoord(block, jj);

          for (site_t kk = 0; kk < GetBlockSize(); kk++)
          {
            site_t site_k = GetSiteCoord(block, kk);

            ++m;

            unsigned int *site_type = &Blocks[block].site_data[m];
            if (!reader->readUnsignedInt(*site_type))
            {
              std::cout << "Error reading site type\n";
            }

            if ( (*site_type & SITE_TYPE_MASK) == SOLID_TYPE)
            {
              Blocks[block].ProcessorRankForEachBlockSite[m] = BIG_NUMBER2;
              continue;
            }

            Blocks[block].ProcessorRankForEachBlockSite[m] = -1;

            if (GetCollisionType(*site_type) != FLUID)
            {
              // Neither solid nor simple fluid
              if (Blocks[block].wall_data == NULL)
              {
                Blocks[block].wall_data = new WallData[GetSitesPerBlockVolumeUnit()];
              }

              if (GetCollisionType(*site_type) & INLET || GetCollisionType(*site_type) & OUTLET)
              {
                double temp;
                // INLET or OUTLET or both.
                // These values are the boundary normal and the boundary distance.
                for (int l = 0; l < 3; l++)
                {
                  if (!reader->readDouble(temp))
                  {
                    std::cout << "Error reading boundary normals\n";
                  }
                }

                if (!reader->readDouble(temp))
                {
                  std::cout << "Error reading boundary distances\n";
                }
              }

              if (GetCollisionType(*site_type) & EDGE)
              {
                // EDGE bit set
                for (int l = 0; l < 3; l++)
                {
                  if (!reader->readDouble(Blocks[block].wall_data[m].wall_nor[l]))
                  {
                    std::cout << "Error reading edge normal\n";
                  }
                }

                double temp;
                if (!reader->readDouble(temp))
                {
                  std::cout << "Error reading edge distance\n";
                }
              }

              for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
              {
                if (!reader->readDouble(Blocks[block].wall_data[m].cut_dist[l]))
                {
                  std::cout << "Error reading cut distances\n";
                }
              }
            }
          } // kk
        } // jj
      } // ii
    }

  }
}
