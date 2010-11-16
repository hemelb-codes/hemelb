#ifndef HEMELB_LB_GLOBALLATTICEDATA_H
#define HEMELB_LB_GLOBALLATTICEDATA_H

#include "D3Q15.h"
#include <cstdlib>

namespace hemelb
{
  namespace lb
  {
    // Data about an element of the domain wall
    struct WallData
    {
        // estimated boundary normal.
        double boundary_nor[3];
        // estimated minimum distance (in lattice units) from the
        // boundary;
        double boundary_dist;
        // estimated wall normal (if the site is close to the wall);
        double wall_nor[3];
        // estimated minimum distance (in lattice units) from the wall;
        // if the site is close to the wall surface
        double wall_dist;
        // cut distances along the 14 non-zero lattice vectors;
        // each one is between 0 and 1 if the surface cuts the corresponding
        // vector or is equal to 1e+30 otherwise
        double cut_dist[D3Q15::NUMVECTORS - 1];
    };

    // Data about each global block in the lattice,
    // site_data[] is an array containing individual lattice site data
    // within a global block.
    struct BlockData
    {
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

    struct GlobalLatticeData
    {
        int SitesX, SitesY, SitesZ;
        int BlocksX, BlocksY, BlocksZ;
        int BlockCount;
        int BlockSize;
        int Log2BlockSize;
        int SitesPerBlockVolumeUnit;
        BlockData * Blocks;

        ~GlobalLatticeData()
        {
          delete[] Blocks;
        }

        // Function that finds the pointer to the rank on which a particular site
        // resides. If the site is in an empty block, return NULL.
        int * GetProcIdFromGlobalCoords(int iSiteI, int iSiteJ, int iSiteK)
        {
          // If the given site location is outside the bounding box return a NULL
          // pointer.
          if (iSiteI < 0 || iSiteI >= SitesX || iSiteJ < 0 || iSiteJ >= SitesY
              || iSiteK < 0 || iSiteK >= SitesZ)
          {
            return NULL;
          }

          // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
          int i = iSiteI >> Log2BlockSize;
          int j = iSiteJ >> Log2BlockSize;
          int k = iSiteK >> Log2BlockSize;

          // Get the block from the block identifiers.
          BlockData * lBlock = &Blocks[ (i * BlocksY + j) * BlocksZ + k];

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

        unsigned int * GetSiteData(int iSiteI, int iSiteJ, int iSiteK)
        {
          // If site is out of the bounding box, return NULL.
          if (iSiteI < 0 || iSiteI >= SitesX || iSiteJ < 0 || iSiteJ >= SitesY
              || iSiteK < 0 || iSiteK >= SitesZ)
          {
            return NULL;
          }

          // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
          int i = iSiteI >> Log2BlockSize;
          int j = iSiteJ >> Log2BlockSize;
          int k = iSiteK >> Log2BlockSize;

          // Pointer to the block
          BlockData * lBlock = &Blocks[ (i * BlocksY + j) * BlocksZ + k];

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
    };
  }
}

#endif /* HEMELB_LB_GLOBALLATTICEDATA_H */
