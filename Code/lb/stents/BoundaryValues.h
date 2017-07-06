
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STENTS_BOUNDARYVALUES_H
#define HEMELB_LB_STENTS_BOUNDARYVALUES_H

#include "net/IOCommunicator.h"
#include "net/IteratedAction.h"
#include "lb/stents/Stent.h"
#include "geometry/LatticeData.h"
#include "lb/iolets/BoundaryCommunicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      class BoundaryValues : public net::IteratedAction
      {
        public:
          BoundaryValues(geometry::LatticeData* latticeData,
                         const std::vector<stents::Stent*> &stents,
                         SimulationState* simulationState,
                         const net::MpiCommunicator& comms,
                         const util::UnitConverter& units);
          ~BoundaryValues();

          void RequestComms();
          void EndIteration();
          void Reset();

          void FinishReceive();

          LatticeDensity GetBoundaryDensity(const int index);

          LatticeDensity GetDensityMin(int boundaryId);
          LatticeDensity GetDensityMax(int boundaryId);

          static proc_t GetBCProcRank();
          stents::Stent* GetLocalStent(unsigned int index)
          {
            return stents[localStentIDs[index]];
          }
          unsigned int GetLocalStentCount()
          {
            return localStentCount;
          }
          inline unsigned int GetTimeStep() const
          {
            return state->GetTimeStep();
          }

        private:
          bool IsStentOnThisProc(geometry::LatticeData* latticeData, int boundaryId);
          std::vector<int> GatherProcList(bool hasBoundary);
          void HandleComms(stents::Stent* stent);
          int totalStentCount;
          // Number of Stents and vector of their indices for communication purposes
          int localStentCount;
          std::vector<int> localStentIDs;
          // Has to be a vector of pointers for InOutLet polymorphism
          std::vector<stents::Stent*> stents;

          SimulationState* state;
          const util::UnitConverter& unitConverter;
          BoundaryCommunicator bcComms;
      }
      ;
    }
  }
}

#endif /* HEMELB_LB_STENTS_BOUNDARYVALUES_H */
