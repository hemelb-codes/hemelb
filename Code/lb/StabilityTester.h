#ifndef HEMELB_LB_STABILITYCHECKER_H
#define HEMELB_LB_STABILITYCHECKER_H

#include "net/PhasedBroadcast.h"
#include "lb/LocalLatticeData.h"

namespace hemelb
{
  namespace lb
  {
    class StabilityTester : public net::PhasedBroadcast
    {
      public:
        StabilityTester(const hemelb::lb::LocalLatticeData * iLocalLatDat, int * bStability);

        virtual void Reset();

      protected:
        virtual void ProgressFromChildren();
        virtual void ProgressFromParent();
        virtual void ProgressToChildren();
        virtual void ProgressToParent();

        virtual void Effect();

      private:
        static const unsigned int SPREADFACTOR = 10;
        const hemelb::lb::LocalLatticeData * mLocalLatDat;

        int mLocalStability;
        int mChildrensStability[SPREADFACTOR];
        int mReceivedSimulationStability;
        int * mPublicSimulationStability;
    };
  }
}

#endif /* HEMELB_LB_STABILITYCHECKER_H */
