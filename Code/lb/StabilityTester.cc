#include "lb/StabilityTester.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {

    StabilityTester::StabilityTester(const geometry::LatticeData * iLatDat,
                                     net::Net* net,
                                     SimulationState* simState,
                                     reporting::Timers& timings) :
        net::PhasedBroadcastRegular<>(net, simState, SPREADFACTOR), mLatDat(iLatDat), mSimState(simState), timings(timings)
    {
      Reset();
    }

    void StabilityTester::Reset()
    {
      // Re-initialise all values to be Stable.
      mUpwardsStability = Stable;
      mDownwardsStability = Stable;

      mSimState->SetStability(Stable);

      for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
      {
        mChildrensStability[ii] = Stable;
      }
    }

    void StabilityTester::PostReceiveFromChildren(unsigned long splayNumber)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();

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

      timings[hemelb::reporting::Timers::monitoring].Stop();
    }

    void StabilityTester::ProgressFromChildren(unsigned long splayNumber)
    {
      ReceiveFromChildren<int>(mChildrensStability, 1);
    }

    void StabilityTester::ProgressFromParent(unsigned long splayNumber)
    {
      ReceiveFromParent<int>(&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToChildren(unsigned long splayNumber)
    {
      SendToChildren<int>(&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToParent(unsigned long splayNumber)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();

      // No need to bother testing out local lattice points if we're going to be
      // sending up a 'Unstable' value anyway.
      if (mUpwardsStability != Unstable)
      {
        for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
        {
          for (unsigned int l = 0; l < mLatDat->GetLatticeInfo().GetNumVectors(); l++)
          {
            distribn_t value = *mLatDat->GetFNew(i * mLatDat->GetLatticeInfo().GetNumVectors() + l);

            // Note that by testing for value > 0.0, we also catch stray NaNs.
            if (! (value > 0.0))
            {
              mUpwardsStability = Unstable;
            }
          }
        }

        timings[hemelb::reporting::Timers::monitoring].Stop();
      }

      SendToParent<int>(&mUpwardsStability, 1);
    }

    void StabilityTester::TopNodeAction()
    {
      mDownwardsStability = mUpwardsStability;
    }

    void StabilityTester::Effect()
    {
      mSimState->SetStability((Stability) mDownwardsStability);
    }

  }
}
