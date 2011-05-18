#include "lb/StabilityTester.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {

    StabilityTester::StabilityTester(const geometry::LatticeData * iLatDat,
                                     net::Net* net,
                                     SimulationState* simState) :
      net::PhasedBroadcastRegular<>(net, simState, SPREADFACTOR), mLatDat(iLatDat)
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

    void StabilityTester::PostReceiveFromChildren(unsigned int splayNumber)
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

    void StabilityTester::ProgressFromChildren(unsigned int splayNumber)
    {
      ReceiveFromChildren<int> (mChildrensStability, 1);
    }

    void StabilityTester::ProgressFromParent(unsigned int splayNumber)
    {
      ReceiveFromParent<int> (&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToChildren(unsigned int splayNumber)
    {
      SendToChildren<int> (&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToParent(unsigned int splayNumber)
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
