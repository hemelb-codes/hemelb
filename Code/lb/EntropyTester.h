// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_ENTROPYTESTER_H
#define HEMELB_LB_ENTROPYTESTER_H

#include "net/CollectiveAction.h"
#include "geometry/LatticeData.h"
#include "lb/HFunction.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace lb
  {
    template<class LatticeType>
    class EntropyTester : public net::CollectiveAction
    {
      public:
        EntropyTester(int* collisionTypes, unsigned int typesTested,
                      const geometry::LatticeData * iLatDat, net::Net* net,
                      SimulationState* simState, reporting::Timers& timings) :
            net::CollectiveAction(net->GetCommunicator(), timings), mLatDat(iLatDat),
                mHPreCollision(mLatDat->GetLocalFluidSiteCount())
        {
          for (unsigned int i = 0; i < COLLISION_TYPES; i++)
          {
            mCollisionTypesTested[i] = false;
          }
          for (unsigned int i = 0; i < typesTested; i++)
          {
            mCollisionTypesTested[collisionTypes[i]] = true;
          }

          Reset();
        }

        // Run at BeginPhase
        void RequestComms(void)
        {
          // Store pre-collision values.
          site_t offset = 0;

          for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
          {
            if (mCollisionTypesTested[collision_type])
            {
              for (site_t i = offset;
                  i < offset + mLatDat->GetMidDomainCollisionCount(collision_type); i++)
              {
                HFunction<LatticeType> HFuncOldOld(mLatDat->GetFNew(LatticeType::NUMVECTORS*i), NULL);
                mHPreCollision[i] = HFuncOldOld.eval();
              }
            }

            offset += mLatDat->GetMidDomainCollisionCount(collision_type);
          }

          for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
          {
            if (mCollisionTypesTested[collision_type])
            {
              for (site_t i = offset;
                  i < offset + mLatDat->GetDomainEdgeCollisionCount(collision_type); i++)
              {
                HFunction<LatticeType> HFuncOldOld(mLatDat->GetFNew(LatticeType::NUMVECTORS*i), NULL);
                mHPreCollision[i] = HFuncOldOld.eval();
              }
            }

            offset += mLatDat->GetDomainEdgeCollisionCount(collision_type);
          }
        }
        void PreSend()
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
              for (site_t i = offset;
                  i < offset + mLatDat->GetMidDomainCollisionCount(collision_type); i++)
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
              for (site_t i = offset;
                  i < offset + mLatDat->GetDomainEdgeCollisionCount(collision_type); i++)
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
            localHTheorem = DISOBEYED;
          }
        }
        void Send(void)
        {
          collectiveReq = collectiveComm.Ireduce(localHTheorem, MPI_MAX, RootRank, globalHTheorem);
        }
        /**
         * Take the combined stability information (an int, with a value of hemelb::lb::Unstable
         * if any child node is unstable) and start passing it back down the tree.
         */
        void PostReceive()
        {
          if (collectiveComm.Rank() == RootRank && globalHTheorem == DISOBEYED)
          {
            log::Logger::Log<log::Error, log::Singleton>("H Theorem violated.");
          }
        }
        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset()
        {
          // Re-initialise all values to indicate that the H-theorem is obeyed.
          localHTheorem = OBEYED;
          globalHTheorem = OBEYED;
        }
      private:
        enum HTHEOREM
        {
          OBEYED = 0,
          DISOBEYED = 1
        };
        static const int RootRank = 0;

        /**
         * Slightly arbitrary spread factor for the tree.
         */
        const geometry::LatticeData * mLatDat;

        int localHTheorem;
        int globalHTheorem;

        bool mCollisionTypesTested[COLLISION_TYPES];
        std::vector<double> mHPreCollision;
    };

  }
}

#endif /* HEMELB_LB_ENTROPYTESTER_H */
