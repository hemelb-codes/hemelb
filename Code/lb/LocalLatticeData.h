#ifndef HEMELB_LB_LOCALLATTICEDATA_H
#define HEMELB_LB_LOCALLATTICEDATA_H

namespace hemelb
{
  namespace lb
  {
    struct LocalLatticeData
    {
        double *FOld;
        double *FNew;
        int *FNeighbours;
    };
  }
}

#endif /* HEMELB_LB_LOCALLATTICEDATA_H */
