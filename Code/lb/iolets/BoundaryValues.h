// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_BOUNDARYVALUES_H
#define HEMELB_LB_IOLETS_BOUNDARYVALUES_H

#include "net/IOCommunicator.h"
#include "net/IteratedAction.h"
#include "lb/iolets/InOutLet.h"
#include "geometry/LatticeData.h"
#include "lb/iolets/BoundaryCommunicator.h"
#include "util/clone_ptr.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class BoundaryValues : public net::IteratedAction
      {
	using IoletPtr = util::clone_ptr<iolets::InOutLet>;
        public:
          BoundaryValues(geometry::SiteType ioletType, geometry::LatticeData* latticeData,
                         const std::vector<IoletPtr>& iolets,
                         SimulationState* simulationState, const net::MpiCommunicator& comms,
                         const util::UnitConverter& units);

          void RequestComms();
          void EndIteration();
          void Reset();

          void FinishReceive();

          LatticeDensity GetBoundaryDensity(const int index);

          LatticeDensity GetDensityMin(int boundaryId);
          LatticeDensity GetDensityMax(int boundaryId);

          static proc_t GetBCProcRank();

	  // Borrow the pointer to an Iolet - this object still owns
	  // the value.
          iolets::InOutLet* GetLocalIolet(unsigned int index)
          {
            return iolets[localIoletIDs[index]].get();
          }
          unsigned int GetLocalIoletCount()
          {
            return localIoletCount;
          }
          inline unsigned int GetTimeStep() const
          {
            return state->GetTimeStep();
          }
          inline geometry::SiteType GetIoletType() const
          {
            return ioletType;
          }

        private:
          bool IsIOletOnThisProc(geometry::SiteType ioletType, geometry::LatticeData* latticeData,
                                 int boundaryId);
          std::vector<int> GatherProcList(bool hasBoundary);
          void HandleComms(iolets::InOutLet* iolet);
          geometry::SiteType ioletType;
          int totalIoletCount;
          // Number of IOlets and vector of their indices for communication purposes
          int localIoletCount;
          std::vector<int> localIoletIDs;
          // Has to be a vector of pointers for InOutLet polymorphism
          std::vector<IoletPtr> iolets;

          SimulationState* state;
          const util::UnitConverter& unitConverter;
          BoundaryCommunicator bcComms;
      }
      ;
    }
  }
}

#endif /* HEMELB_LB_IOLETS_BOUNDARYVALUES_H */
