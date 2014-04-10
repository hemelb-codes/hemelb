// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_IOLETS_BOUNDARYVALUES_H
#define HEMELB_LB_IOLETS_BOUNDARYVALUES_H

#include "net/IOCommunicator.h"
#include "net/IteratedAction.h"
#include "lb/iolets/InOutLet.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class BoundaryValues : public net::IteratedAction
      {
        public:
          BoundaryValues(geometry::SiteType ioletType,
                         geometry::LatticeData* latticeData,
                         const std::vector<iolets::InOutLet*> &iolets,
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

          static bool IsCurrentProcTheBCProc();
          static proc_t GetBCProcRank();
          iolets::InOutLet* GetLocalIolet(unsigned int index)
          {
            return iolets[localIoletIDs[index]];
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
          bool IsIOletOnThisProc(geometry::SiteType ioletType, geometry::LatticeData* latticeData, int boundaryId);
          std::vector<int> GatherProcList(bool hasBoundary);
          void HandleComms(iolets::InOutLet* iolet);
          geometry::SiteType ioletType;
          int totalIoletCount;
          // Number of IOlets and vector of their indices for communication purposes
          int localIoletCount;
          std::vector<int> localIoletIDs;
          // Has to be a vector of pointers for InOutLet polymorphism
          std::vector<iolets::InOutLet*> iolets;

          SimulationState* state;
          const util::UnitConverter& unitConverter;
          net::MpiCommunicator privateComms;
      }
      ;
    }
  }
}

#endif /* HEMELB_LB_IOLETS_BOUNDARYVALUES_H */
