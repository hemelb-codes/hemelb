#ifndef HEMELB_LB_BOUNDARIES_BOUNDARYVALUES_H
#define HEMELB_LB_BOUNDARIES_BOUNDARYVALUES_H

#include "lb/boundaries/BoundaryComms.h"
#include "topology/NetworkTopology.h"
#include "net/IteratedAction.h"
#include "lb/boundaries/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      class BoundaryValues : public net::IteratedAction
      {
        public:
          BoundaryValues(geometry::LatticeData::SiteType IOtype,
                         geometry::LatticeData* iLatDat,
                         std::vector<iolets::InOutLet*> &iiolets,
                         SimulationState* iSimState,
                         util::UnitConverter* units);
          ~BoundaryValues();

          void RequestComms();
          void EndIteration();
          void Reset();

          void FinishReceive();

          distribn_t GetBoundaryDensity(const int index);

          distribn_t GetDensityMin(int iBoundaryId);
          distribn_t GetDensityMax(int iBoundaryId);

          static bool IsCurrentProcTheBCProc();
          static proc_t GetBCProcRank();

        private:
          std::vector<BoundaryComms*> mComms;
          bool IsIOletOnThisProc(geometry::LatticeData::SiteType IOtype,
                                 geometry::LatticeData* iLatDat,
                                 int iBoundaryId);
          std::vector<int> GatherProcList(bool hasBoundary);

          int nTotIOlets;
          // Number of IOlets and vector of their indices for communication purposes
          int nIOlets;
          std::vector<int> ioletIDs;
          // Has to be a vector of pointers for InOutLet polymorphism
          std::vector<iolets::InOutLet*> iolets;

          std::vector<distribn_t>* density_cycle;

          SimulationState* mState;
          util::UnitConverter* mUnits;
      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_BOUNDARYVALUES_H */
