#include "geometry/LatticeData.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace geometry
  {
    void LatticeData::GlobalLatticeData::SetBasicDetails(site_t iBlocksX,
                                                         site_t iBlocksY,
                                                         site_t iBlocksZ,
                                                         site_t iBlockSize,
                                                         distribn_t iVoxelSize,
                                                         distribn_t iOriginX,
                                                         distribn_t iOriginY,
                                                         distribn_t iOriginZ)
    {
      mBlocksX = iBlocksX;
      mBlocksY = iBlocksY;
      mBlocksZ = iBlocksZ;
      mBlockSize = iBlockSize;
      mVoxelSize = iVoxelSize;
      mOriginX = iOriginX;
      mOriginY = iOriginY;
      mOriginZ = iOriginZ;

      mSitesX = mBlocksX * mBlockSize;
      mSitesY = mBlocksY * mBlockSize;
      mSitesZ = mBlocksZ * mBlockSize;

      mSitesPerBlockVolumeUnit = mBlockSize * mBlockSize * mBlockSize;

      // A shift value we'll need later = log_2(block_size)
      site_t i = mBlockSize;
      log2BlockSize = 0;
      while (i > 1)
      {
        i >>= 1;
        ++log2BlockSize;
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
    distribn_t LatticeData::GlobalLatticeData::GetVoxelSize() const
    {
      return mVoxelSize;
    }
    distribn_t LatticeData::GlobalLatticeData::GetXOrigin() const
    {
      return mOriginX;
    }
    distribn_t LatticeData::GlobalLatticeData::GetYOrigin() const
    {
      return mOriginY;
    }
    distribn_t LatticeData::GlobalLatticeData::GetZOrigin() const
    {
      return mOriginZ;
    }

    site_t LatticeData::GlobalLatticeData::GetSitesPerBlockVolumeUnit() const
    {
      return mSitesPerBlockVolumeUnit;
    }

    bool LatticeData::GlobalLatticeData::IsValidBlock(site_t i, site_t j, site_t k) const
    {
      return i < mBlocksX && j < mBlocksY && k < mBlocksZ && i > -1 && j > -1 && k > -1;
    }

    bool LatticeData::GlobalLatticeData::IsValidLatticeSite(site_t i, site_t j, site_t k) const
    {
      return i >= 0 && j >= 0 && k >= 0 && i < mSitesX && j < mSitesY && k < mSitesZ;
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
    const proc_t* LatticeData::GlobalLatticeData::GetProcIdFromGlobalCoords(const util::Vector3D<
        site_t>& globalSiteCoords) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      site_t blockI = globalSiteCoords.x >> log2BlockSize;
      site_t blockJ = globalSiteCoords.y >> log2BlockSize;
      site_t blockK = globalSiteCoords.z >> log2BlockSize;

      // Get the block from the block identifiers.
      BlockData * block = &Blocks[GetBlockIdFromBlockCoords(blockI, blockJ, blockK)];

      // If an empty (solid) block is addressed, return a NULL pointer.
      if (block->ProcessorRankForEachBlockSite == NULL)
      {
        return NULL;
      }
      else
      {
        // Find site coordinates within the block
        site_t siteI = globalSiteCoords.x - (blockI << log2BlockSize);
        site_t siteJ = globalSiteCoords.y - (blockJ << log2BlockSize);
        site_t siteK = globalSiteCoords.z - (blockK << log2BlockSize);

        // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
        // mProcessorsForEachBlock)
        return &block->ProcessorRankForEachBlockSite[ ( ( (siteI << log2BlockSize) + siteJ)
            << log2BlockSize) + siteK];
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

    void LatticeData::GlobalLatticeData::GetBlockIJK(site_t block,
                                                     site_t* blockI,
                                                     site_t* blockJ,
                                                     site_t* blockK) const
    {
      *blockK = block % GetZBlockCount();
      site_t blockIJData = block / GetZBlockCount();
      *blockJ = blockIJData % GetYBlockCount();
      *blockI = blockIJData / GetYBlockCount();
    }

    site_t LatticeData::GlobalLatticeData::GetSiteCoord(site_t block, site_t localSiteCoord) const
    {
      return (block << log2BlockSize) + localSiteCoord;
    }

    const util::Vector3D<site_t> LatticeData::GlobalLatticeData::GetGlobalCoords(site_t blockNumber,
                                                                                 const util::Vector3D<
                                                                                     site_t>& localSiteCoords) const
    {
      site_t blockI, blockJ, blockK;
      GetBlockIJK(blockNumber, &blockI, &blockJ, &blockK);

      return util::Vector3D<site_t>(GetSiteCoord(blockI, localSiteCoords.x),
                                    GetSiteCoord(blockJ, localSiteCoords.y),
                                    GetSiteCoord(blockK, localSiteCoords.z));
    }

    unsigned int LatticeData::GlobalLatticeData::GetSiteData(site_t siteI,
                                                             site_t siteJ,
                                                             site_t siteK) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      site_t blockI = siteI >> log2BlockSize;
      site_t blockJ = siteJ >> log2BlockSize;
      site_t blockK = siteK >> log2BlockSize;

      // Pointer to the block
      BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(blockI, blockJ, blockK)];

      // Find site coordinates within the block
      site_t localSiteI = siteI - (blockI << log2BlockSize);
      site_t localSiteJ = siteJ - (blockJ << log2BlockSize);
      site_t localSiteK = siteK - (blockK << log2BlockSize);

      // Return pointer to site_data[site]
      return lBlock->site_data[ ( ( (localSiteI << log2BlockSize) + localSiteJ) << log2BlockSize)
          + localSiteK];
    }

    void LatticeData::GlobalLatticeData::ReadBlock(const site_t block, io::writers::xdr::XdrReader* reader)
    {
      if (Blocks[block].site_data == NULL)
      {
        Blocks[block].site_data = new unsigned int[GetSitesPerBlockVolumeUnit()];
      }
      if (Blocks[block].ProcessorRankForEachBlockSite == NULL)
      {
        Blocks[block].ProcessorRankForEachBlockSite = new proc_t[GetSitesPerBlockVolumeUnit()];
      }

      site_t localSiteIndex = -1;
      site_t blockI, blockJ, blockK;

      GetBlockIJK(block, &blockI, &blockJ, &blockK);

      for (site_t localSiteI = 0; localSiteI < GetBlockSize(); localSiteI++)
      {
        site_t globalSiteI = GetSiteCoord(block, localSiteI);

        for (site_t localSiteJ = 0; localSiteJ < GetBlockSize(); localSiteJ++)
        {
          site_t globalSiteJ = GetSiteCoord(block, localSiteJ);

          for (site_t localSiteK = 0; localSiteK < GetBlockSize(); localSiteK++)
          {
            site_t globalSiteK = GetSiteCoord(block, localSiteK);

            ++localSiteIndex;

            unsigned int *site_type = &Blocks[block].site_data[localSiteIndex];
            if (!reader->readUnsignedInt(*site_type))
            {
              std::cout << "Error reading site type\n";
            }

            if ( (*site_type & SITE_TYPE_MASK) == SOLID_TYPE)
            {
              Blocks[block].ProcessorRankForEachBlockSite[localSiteIndex] = BIG_NUMBER2;
              continue;
            }

            Blocks[block].ProcessorRankForEachBlockSite[localSiteIndex] = -1;

            if (GetCollisionType(*site_type) != FLUID)
            {
              // Neither solid nor simple fluid
              if (Blocks[block].wall_data == NULL)
              {
                Blocks[block].wall_data = new WallData[GetSitesPerBlockVolumeUnit()];

                // Initialise all the WallData objects to indicate a lack of walls / cuts.
                for (unsigned int localSiteCount = 0; localSiteCount < GetSitesPerBlockVolumeUnit(); ++localSiteCount)
                {
                  for (unsigned int cutDirection = 0; cutDirection < D3Q15::NUMVECTORS - 1; ++cutDirection)
                  {
                    Blocks[block].wall_data[localSiteCount].cut_dist[cutDirection] = -1.0F;
                  }
                  for (unsigned int wallNormalDimension = 0; wallNormalDimension < 3; ++wallNormalDimension)
                  {
                    Blocks[block].wall_data[localSiteCount].wall_nor[wallNormalDimension] = -1.0F;
                  }
                }
              }

              if (GetCollisionType(*site_type) & INLET || GetCollisionType(*site_type) & OUTLET)
              {
                double temp;
                // INLET or OUTLET or both.
                // These values are the boundary normal and the boundary distance.
                for (int dimension = 0; dimension < 3; dimension++)
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
                for (int dimension = 0; dimension < 3; dimension++)
                {
                  if (!reader->readDouble(Blocks[block].wall_data[localSiteIndex].wall_nor[dimension]))
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

              for (unsigned int direction = 0; direction < (D3Q15::NUMVECTORS - 1); direction++)
              {
                if (!reader->readDouble(Blocks[block].wall_data[localSiteIndex].cut_dist[direction]))
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
