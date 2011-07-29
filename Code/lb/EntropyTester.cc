#include <cstdio>

#include "log/Logger.h"
#include "lb/EntropyTester.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "lb/collisions/implementations/HFunction.h"
using hemelb::lb::collisions::implementations::HFunction;

namespace hemelb
{
  namespace lb
  {

    EntropyTester::EntropyTester(int* collisionTypes,
                                 unsigned int typesTested,
                                 const geometry::LatticeData * iLatDat,
                                 net::Net* net,
                                 SimulationState* simState) :
      net::PhasedBroadcastRegular<false, 1, 1, false, true>(net, simState, SPREADFACTOR),
          mLatDat(iLatDat), mSimState(simState)
    {
      mCollisionTypeTested = new bool[COLLISION_TYPES];
      for (unsigned int i = 0; i < COLLISION_TYPES; i++)
      {
        mCollisionTypeTested[i] = false;
      }
      for (unsigned int i = 0; i < typesTested; i++)
      {
        mCollisionTypeTested[collisionTypes[i]] = true;
      }
      mHPreCollision = new double[mLatDat->GetLocalFluidSiteCount()];

      Reset();
    }

    EntropyTester::~EntropyTester()
    {
      delete[] mCollisionTypeTested;
      delete[] mHPreCollision;
    }

    void EntropyTester::Reset()
    {
      // Re-initialise all values to be Stable.
      mUpwardsStability = Stable;

      mSimState->SetStability(Stable);

      for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
      {
        mChildrensStability[ii] = Stable;
      }
    }

    void EntropyTester::PreReceive()
    {
      // No need to bother testing out local lattice points if we're going to be
      // sending up a 'Unstable' value anyway.
      if (mUpwardsStability != Unstable)
      {
        site_t offset = 0;
        double dH = 0.0;

        // The order of arguments in max is important
        // If distributions go negative HFunc.eval() will return a NaN. Due to the
        // nature of NaN and the structure of max, dH will be assigned as NaN if this is the case
        // This is what we want, because the EntropyTester will fail otherwise and abort when
        // it is simply sufficient to wait until StabilityTester restarts.

        for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
        {
          if (mCollisionTypeTested[collision_type])
          {
            for (site_t i = offset; i < offset + mLatDat->GetInnerCollisionCount(collision_type); i++)
            {
              HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
              dH = util::NumericalFunctions::max(dH, HFunc.eval() - mHPreCollision[i]);
            }
          }

          offset += mLatDat->GetInnerCollisionCount(collision_type);
        }

        for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
        {
          if (mCollisionTypeTested[collision_type])
          {
            for (site_t i = offset; i < offset + mLatDat->GetInterCollisionCount(collision_type); i++)
            {
              HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
              dH = util::NumericalFunctions::max(dH, HFunc.eval() - mHPreCollision[i]);
            }
          }

          offset += mLatDat->GetInterCollisionCount(collision_type);
        }

        /*
         * Ideally dH should never be greater than zero. However, because H is positive definite accuracy as well as
         * rounding and truncation errors can make dH greater than zero in certain cases. The tolerance
         * is limited by the accuracy of the simulation (including the accuracy to which alpha is calculated)
         */
        if (dH > 1.0E-3)
        {
          mUpwardsStability = Unstable;
        }
      }
    }

    void EntropyTester::PostReceiveFromChildren(unsigned long splayNumber)
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

    void EntropyTester::ProgressFromChildren(unsigned long splayNumber)
    {
      ReceiveFromChildren<int> (mChildrensStability, 1);
    }

    void EntropyTester::ProgressToParent(unsigned long splayNumber)
    {
      // Store pre-collision values. Don't bother if unstable already
      if (mUpwardsStability != Unstable)
      {
        site_t offset = 0;

        for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
        {
          if (mCollisionTypeTested[collision_type])
          {
            for (site_t i = offset; i < offset + mLatDat->GetInnerCollisionCount(collision_type); i++)
            {
              HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
              mHPreCollision[i] = HFunc.eval();
            }
          }

          offset += mLatDat->GetInnerCollisionCount(collision_type);
        }

        for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
        {
          if (mCollisionTypeTested[collision_type])
          {
            for (site_t i = offset; i < offset + mLatDat->GetInterCollisionCount(collision_type); i++)
            {
              HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
              mHPreCollision[i] = HFunc.eval();
            }
          }

          offset += mLatDat->GetInterCollisionCount(collision_type);
        }
      }

      SendToParent<int> (&mUpwardsStability, 1);
    }

    void EntropyTester::TopNodeAction()
    {
      if (mUpwardsStability == Unstable)
      {
        std::cout << "!H Theorem violated!" << std::endl;
        int err = MPI_Abort(MPI_COMM_WORLD, 1);

        // This gives us something to work from when we have an error - we get the rank
        // that calls abort, and we get a stack-trace from the exception having been thrown.
        hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Aborting");
        exit(1);
      }
    }

  }
}
