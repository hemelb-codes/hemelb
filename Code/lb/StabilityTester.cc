#include "lb/StabilityTester.h"

namespace hemelb
{
  namespace lb
  {

    StabilityTester::StabilityTester(const hemelb::lb::LocalLatticeData * iLocalLatDat,
                                     int * bStability)
    {
      PhasedBroadcast(SPREADFACTOR);

      mLocalLatDat = iLocalLatDat;
      mPublicSimulationStability = bStability;

      Reset();
    }

    void StabilityTester::Reset()
    {
      mLocalStability = Stable;
      mReceivedSimulationStability = Stable;
      *mPublicSimulationStability = Stable;
      for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
      {
        mChildrensStability[ii] = Stable;
      }
    }

    void StabilityTester::ProgressFromChildren()
    {
      ReceiveFromChildren<int> (mChildrensStability, 1);
    }

    void StabilityTester::ProgressFromParent()
    {
      ReceiveFromParent<int> (&mReceivedSimulationStability, 1);
    }

    void StabilityTester::ProgressToChildren()
    {
      SendToChildren<int> (&mReceivedSimulationStability, 1);
    }

    void StabilityTester::ProgressToParent()
    {
      // No need to test children's stability if this node is already unstable.
      if (mLocalStability == Stable)
      {
        for (int ii = 0; ii < mChildrensStability; ii++)
        {
          if (mChildrensStability == Unstable)
          {
            mLocalStability = Unstable;
            break;
          }
        }
      }

      SendToParent<int> (&mLocalStability, 1);
    }

    void StabilityTester::Effect()
    {
      *mPublicSimulationStability = mReceivedSimulationStability;
    }

  }
