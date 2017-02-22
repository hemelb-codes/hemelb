
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_BOUNDARYVALUES_H
#define HEMELB_LB_IOLETS_BOUNDARYVALUES_H

#include "comm/Communicator.h"
#include "timestep/Actor.h"
#include "lb/iolets/InOutLet.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class BoundaryValues : public timestep::Actor
      {
        public:
          BoundaryValues(geometry::SiteType ioletType,
                         geometry::LatticeData* latticeData,
                         const std::vector<iolets::InOutLet*> &iolets,
                         SimulationState* simulationState,
                         comm::Async::Ptr commQ,
                         const util::UnitConverter& units);
          ~BoundaryValues();

	// timestep::Actor interface
	virtual void BeginAll() { /* not needed */ }
	virtual void Begin();
	virtual void Receive();
	virtual void PreSend() { /* not needed */ }
	virtual void Send();
	virtual void PreWait() { /* not needed */ }
	virtual void Wait() { /* not needed */ }
	virtual void End();
	virtual void EndAll() { /* not needed */ }


	void ForceCommunication();
          // void RequestComms();
          // void EndIteration();
          void Reset();

          // void FinishReceive();

          LatticeDensity GetBoundaryDensity(const int index);

  	  static constexpr int GetBCProcRank() {
	    return 0;
	  }

	  inline bool IsCurrentProcTheBCProc() {
	    return asyncCommsQ->GetComm()->Rank() == GetBCProcRank();
	  }

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

	inline const std::vector<int>& GetProcsForIolet(int i) const {
	  return procsForIolet[i];
	}
	inline const std::vector<int>& GetProcsForIolet(const InOutLet* io) const {
	  auto iter = std::find(iolets.begin(), iolets.end(), io);
	  if (iter == iolets.end())
	    throw Exception() << "Can't find iolet in BoundaryValues";
	  return  GetProcsForIolet(iter - iolets.begin());
	}
	
        private:
          bool IsIOletOnThisProc(geometry::SiteType ioletType, geometry::LatticeData* latticeData, int boundaryId);
          std::vector<int> GatherProcList(bool hasBoundary);
          // void HandleComms(iolets::InOutLet* iolet);
          geometry::SiteType ioletType;
          int totalIoletCount;
          // Number of IOlets and vector of their indices for communication purposes
          int localIoletCount;
          std::vector<int> localIoletIDs;
          // Has to be a vector of pointers for InOutLet polymorphism
          std::vector<iolets::InOutLet*> iolets;
	  // For each iolet, lists the procs that have it (only populated on BCProc rank)
	  std::vector<std::vector<int>> procsForIolet;
          SimulationState* state;
          const util::UnitConverter& unitConverter;
          comm::Async::Ptr asyncCommsQ;
      }
      ;
    }
  }
}

#endif /* HEMELB_LB_IOLETS_BOUNDARYVALUES_H */
