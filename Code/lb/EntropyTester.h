
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_ENTROPYTESTER_H
#define HEMELB_LB_ENTROPYTESTER_H

#include "net/PhasedBroadcastRegular.h"
#include "geometry/LatticeData.h"
#include "lb/HFunction.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace lb
  {
    template<class LatticeType>
    class EntropyTester : public net::PhasedBroadcastRegular<false, 1, 1, false, true>
    {
      public:
        EntropyTester(int* collisionTypes,
                      unsigned int typesTested,
                      const geometry::LatticeData * iLatDat,
                      net::Net* net,
                      SimulationState* simState) :
            net::PhasedBroadcastRegular<false, 1, 1, false, true>(net, simState, SPREADFACTOR), mLatDat(iLatDat)
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

        ~EntropyTester()
        {
          delete[] mHPreCollision;
        }

        void PreReceive()
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
              for (site_t i = offset; i < offset + mLatDat->GetMidDomainCollisionCount(collision_type); i++)
              {
                const geometry::Site<const geometry::LatticeData> site = mLatDat->GetSite(i);

                HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), NULL);
                dHMax = util::NumericalFunctions::max(dHMax, HFunc.eval() - mHPreCollision[i]);
              }
            }

            offset += mLatDat->GetMidDomainCollisionCount(collision_type);
          }

          for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
          {
            if (mCollisionTypesTested[collision_type])
            {
              for (site_t i = offset; i < offset + mLatDat->GetDomainEdgeCollisionCount(collision_type); i++)
              {
                const geometry::Site<const geometry::LatticeData> site = mLatDat->GetSite(i);

                HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), NULL);
                dHMax = util::NumericalFunctions::max(dHMax, HFunc.eval() - mHPreCollision[i]);
              }
            }

            offset += mLatDat->GetDomainEdgeCollisionCount(collision_type);
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

        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset()
        {
          // Re-initialise all values to indicate that the H-theorem is obeyed.
          mUpwardsValue = OBEYED;

          for (unsigned int ii = 0; ii < SPREADFACTOR; ii++)
          {
            mChildrensValues[ii] = OBEYED;
          }
        }

      protected:
        /**
         * Override the methods from the base class to propagate data from the root, and
         * to send data about this node and its childrens' stabilities up towards the root.
         */
        void ProgressFromChildren(unsigned long splayNumber)
        {
          ReceiveFromChildren<int>(mChildrensValues, 1);
        }

        void ProgressToParent(unsigned long splayNumber)
        {
          // Store pre-collision values.
          site_t offset = 0;

          for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
          {
            if (mCollisionTypesTested[collision_type])
            {
              for (site_t i = offset; i < offset + mLatDat->GetMidDomainCollisionCount(collision_type); i++)
              {
                const geometry::Site<const geometry::LatticeData> site = mLatDat->GetSite(i);
                HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), NULL);
                mHPreCollision[i] = HFunc.eval();
              }
            }

            offset += mLatDat->GetMidDomainCollisionCount(collision_type);
          }

          for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
          {
            if (mCollisionTypesTested[collision_type])
            {
              for (site_t i = offset; i < offset + mLatDat->GetDomainEdgeCollisionCount(collision_type); i++)
              {
                const geometry::Site<const geometry::LatticeData> site = mLatDat->GetSite(i);
                HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), NULL);
                mHPreCollision[i] = HFunc.eval();
              }
            }

            offset += mLatDat->GetDomainEdgeCollisionCount(collision_type);
          }

          SendToParent<int>(&mUpwardsValue, 1);
        }

        /**
         * Take the combined stability information (an int, with a value of hemelb::lb::Unstable
         * if any child node is unstable) and start passing it back down the tree.
         */
        void TopNodeAction()
        {
          if (mUpwardsValue == DISOBEYED)
          {
            log::Logger::Log<log::Error, log::Singleton>("H Theorem violated.");
          }
        }

        /**
         * Override the method from the base class to use the data from child nodes.
         */
        void PostReceiveFromChildren(unsigned long splayNumber)
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

      private:
        enum HTHEOREM
        {
          OBEYED,
          DISOBEYED
        };

        /**
         * Slightly arbitrary spread factor for the tree.
         */
        static const unsigned int SPREADFACTOR = 10;

        const geometry::LatticeData * mLatDat;

        /**
         * Stability value of this node and its children to propagate upwards.
         */
        int mUpwardsValue;
        /**
         * Array for storing the passed-up stability values from child nodes.
         */
        int mChildrensValues[SPREADFACTOR];

        bool mCollisionTypesTested[COLLISION_TYPES];
        double* mHPreCollision;
    };

  }
}

#endif /* HEMELB_LB_ENTROPYTESTER_H */
