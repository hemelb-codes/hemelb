#include "lb/StabilityTester.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {

    StabilityTester::StabilityTester(const geometry::LatticeData * iLatDat,
                                     net::Net* net,
                                     topology::NetworkTopology* iNetTop,
                                     SimulationState* simState) :
      net::PhasedBroadcast(net, iNetTop, simState, SPREADFACTOR), mLatDat(iLatDat)
    {
      mPublicSimulationStability = &simState->Stability;
      Reset();
    }

    void StabilityTester::Reset()
    {
      // Re-initialise all values to be Stable.
      mUpwardsStability = Stable;
      mDownwardsStability = Stable;
      *mPublicSimulationStability = Stable;
      for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
      {
        mChildrensStability[ii] = Stable;
      }
    }

    void StabilityTester::PostReceiveFromChildren()
    {
      // No need to test children's stability if this node is already unstable.
      if (mUpwardsStability == Stable)
      {
        for (int ii = 0; ii < (int) SPREADFACTOR; ii++)
        {
          if (mChildrensStability[ii] == Unstable)
          {
            mUpwardsStability = Unstable;
            break;
          }
        }
      }
    }

    void StabilityTester::ProgressFromChildren()
    {
      ReceiveFromChildren<int> (mChildrensStability, 1);
    }

    void StabilityTester::ProgressFromParent()
    {
      ReceiveFromParent<int> (&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToChildren()
    {
      SendToChildren<int> (&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToParent()
    {
      // No need to bother testing out local lattice points if we're going to be
      // sending up a 'Unstable' value anyway.
      if (mUpwardsStability != Unstable)
      {
        for (unsigned int i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
        {
          for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
          {
            if (*mLatDat->GetFNew(i * D3Q15::NUMVECTORS + l) < 0.0)
            {
              mUpwardsStability = Unstable;
            }
          }
        }
      }

      SendToParent<int> (&mUpwardsStability, 1);
    }

    void StabilityTester::TopNodeAction()
    {
      mDownwardsStability = mUpwardsStability;
    }

    void StabilityTester::Effect()
    {
      *mPublicSimulationStability = mDownwardsStability;
    }

  }
}
