#include <cstdio>

#include "lb/EntropyTester.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "lb/collisions/implementations/HFunction.h"
using hemelb::lb::collisions::implementations::HFunction;

namespace hemelb
{
  namespace lb
  {

    EntropyTester::EntropyTester(unsigned int* collisionTypes,
                                 unsigned int typesTested,
                                 const geometry::LatticeData * iLatDat,
                                 net::Net* net,
                                 SimulationState* simState) :
      net::PhasedBroadcastRegular<>(net, simState, SPREADFACTOR), mLatDat(iLatDat),
          mSimState(simState)
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
      mDownwardsStability = Stable;

      mSimState->SetStability(Stable);

      for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
      {
        mChildrensStability[ii] = Stable;
      }
    }

    void EntropyTester::RequestComms()
    {
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

      const unsigned long iCycleNumber = Get0IndexedIterationNumber();
      const unsigned long firstAscent = base::GetFirstAscending();
      const unsigned long firstDescent = base::GetFirstDescending();
      const unsigned long traversalLength = base::GetTraverseTime();

      // Nothing to do for initial action case.

      // Next, deal with the case of a cycle with an initial pass down the tree.
      if (iCycleNumber >= firstDescent && iCycleNumber < firstAscent)
      {
        unsigned long sendOverlap;
        unsigned long receiveOverlap;

        if (base::GetSendChildrenOverlap(iCycleNumber - firstDescent, &sendOverlap))
        {
          ProgressToChildren(sendOverlap);
        }

        if (base::GetReceiveParentOverlap(iCycleNumber - firstDescent, &receiveOverlap))
        {
          ProgressFromParent(receiveOverlap);
        }
      }

      // And deal with the case of a cycle with a pass up the tree.
      if (iCycleNumber >= firstAscent)
      {
        unsigned long sendOverlap;
        unsigned long receiveOverlap;

        if (base::GetSendParentOverlap(iCycleNumber - firstAscent, &sendOverlap))
        {
          ProgressToParent(sendOverlap);
        }

        if (base::GetReceiveChildrenOverlap(iCycleNumber - firstAscent, &receiveOverlap))
        {
          ProgressFromChildren(receiveOverlap);
        }
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

        for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
        {
          if (mCollisionTypeTested[collision_type])
          {
            for (site_t i = offset; i < offset + mLatDat->GetInnerCollisionCount(collision_type); i++)
            {
              HFunction HFunc(mLatDat->GetFOld(i * D3Q15::NUMVECTORS), NULL);
              dH = util::NumericalFunctions::max(HFunc.eval() - mHPreCollision[i], dH);
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
              dH = util::NumericalFunctions::max(HFunc.eval() - mHPreCollision[i], dH);
            }
          }

          offset += mLatDat->GetInterCollisionCount(collision_type);
        }

        if (dH > 1.0E-10)
        {
          mUpwardsStability = Unstable;
          // If mDownwardsStability is also set to unstable program freezes...
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

    void EntropyTester::ProgressFromParent(unsigned long splayNumber)
    {
      ReceiveFromParent<int> (&mDownwardsStability, 1);
    }

    void EntropyTester::ProgressToChildren(unsigned long splayNumber)
    {
      SendToChildren<int> (&mDownwardsStability, 1);
    }

    void EntropyTester::ProgressToParent(unsigned long splayNumber)
    {
      SendToParent<int> (&mUpwardsStability, 1);
    }

    void EntropyTester::TopNodeAction()
    {
      mDownwardsStability = mUpwardsStability;
    }

    void EntropyTester::Effect()
    {
      mSimState->SetStability((Stability) mDownwardsStability);
    }

  }
}
