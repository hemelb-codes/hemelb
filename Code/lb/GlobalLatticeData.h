#ifndef HEMELB_LB_GLOBALLATTICEDATA_H
#define HEMELB_LB_GLOBALLATTICEDATA_H

namespace hemelb
{
  namespace lb
  {
    struct GlobalLatticeData
    {
        int SitesX, SitesY, SitesZ;
        int BlocksX, BlocksY, BlocksZ;
        int BlockCount;
        int BlockSize;
        int Log2BlockSize;
        int SitesPerBlockVolumeUnit;
    };
  }
}

#endif /* HEMELB_LB_GLOBALLATTICEDATA_H */
