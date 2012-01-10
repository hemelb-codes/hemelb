#include <cstdio>

#include "log/Logger.h"
#include "lb/EntropyTester.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "lb/HFunction.h"

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
          mLatDat(iLatDat)
    {
      for (unsigned int i = 0; i < COLLISION_TYPES; i++)
      {
        mCollisionTypesTested[i] = false;
      }
      for (unsigned int i = 0; i < typesTested; i++)
      {
        mCollisionTypesTested[collisionTypes[i]] = true;
      }

      mHPreCollision = new double[mLatDat->GetLocalFluidSiteCount()];

      Reset();
    }

    EntropyTester::~EntropyTester()
    {
      delete[] mHPreCollision;
    }

    void EntropyTester::Reset()
    {
      // Re-initialise all values to indicate that the H-theorem is obeyed.
      mUpwardsValue = OBEYED;

      for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
      {
        mChildrensValues[ii] = OBEYED;
      }
    }

    void EntropyTester::PreReceive()
    {
      site_t offset = 0;
      double dHMax = 0.0;

      // The order of arguments in max is important
      // If distributions go negative HFunc.eval() will return a NaN. Due to the
      // nature of NaN and the structure of max, dH will be assigned as NaN if this is the case
      // This is what we want, because the EntropyTester will fail otherwise and abort when
      // it is simply sufficient to wait until StabilityTester restarts.

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        if (mCollisionTypesTested[collision_type])
        {
          for (site_t i = offset; i < offset + mLatDat->GetInnerCollisionCount(collision_type); i++)
          {
            HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
            dHMax = util::NumericalFunctions::max(dHMax, HFunc.eval() - mHPreCollision[i]);
          }
        }

        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        if (mCollisionTypesTested[collision_type])
        {
          for (site_t i = offset; i < offset + mLatDat->GetInterCollisionCount(collision_type); i++)
          {
            HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
            dHMax = util::NumericalFunctions::max(dHMax, HFunc.eval() - mHPreCollision[i]);
          }
        }

        offset += mLatDat->GetInterCollisionCount(collision_type);
      }

      /*
       * Ideally dH should never be greater than zero. However, because H is positive definite accuracy as well as
       * rounding and truncation errors can make dH greater than zero in certain cases. The tolerance
       * is limited by the accuracy of the simulation (including the accuracy to which alpha is calculated)
       * The tolerance has to be at least as big as the accuracy to which alpha is calculated
       */
      if (dHMax > 1.0E-6)
      {
        mUpwardsValue = DISOBEYED;
      }
    }

    void EntropyTester::PostReceiveFromChildren(unsigned long splayNumber)
    {
      // No need to test children's entropy direction if this node already disobeys H-theorem.
      if (mUpwardsValue == OBEYED)
      {
        for (int ii = 0; ii < (int) SPREADFACTOR; ii++)
        {
          if (mChildrensValues[ii] == DISOBEYED)
          {
            mUpwardsValue = DISOBEYED;
            break;
          }
        }
      }
    }

    void EntropyTester::ProgressFromChildren(unsigned long splayNumber)
    {
      ReceiveFromChildren<int> (mChildrensValues, 1);
    }

    void EntropyTester::ProgressToParent(unsigned long splayNumber)
    {
      // Store pre-collision values.
      site_t offset = 0;

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        if (mCollisionTypesTested[collision_type])
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
        if (mCollisionTypesTested[collision_type])
        {
          for (site_t i = offset; i < offset + mLatDat->GetInterCollisionCount(collision_type); i++)
          {
            HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
            mHPreCollision[i] = HFunc.eval();
          }
        }

        offset += mLatDat->GetInterCollisionCount(collision_type);
      }

      SendToParent<int> (&mUpwardsValue, 1);
    }

    void EntropyTester::TopNodeAction()
    {
      if (mUpwardsValue == DISOBEYED)
      {
        log::Logger::Log<log::Info, log::Singleton>("H Theorem violated.");
      }
    }

  }
}
