#ifndef HEMELB_LB_BOUNDARIES_BOUNDARYVALUES_H
#define HEMELB_LB_BOUNDARIES_BOUNDARYVALUES_H

#include "lb/boundaries/BoundaryComms.h"
#include "topology/NetworkTopology.h"
#include "net/IteratedAction.h"

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
                         SimConfig* iSimConfig,
                         SimulationState* iSimState,
                         util::UnitConverter* iUnits);
          ~BoundaryValues();

          void RequestComms();
          void EndIteration();
          void Reset();

          void FinishReceive();

          void ResetPrePeriodChange();
          void ResetPostPeriodChange();

          distribn_t GetBoundaryDensity(const int index);

          distribn_t GetDensityMin(int iBoundaryId);
          distribn_t GetDensityMax(int iBoundaryId);

          static bool IsCurrentProcTheBCProc();

        private:
          proc_t BCproc;
          void FindBCProcRank();

          BoundaryComms** mComms;
          bool IsIOletOnThisProc(geometry::LatticeData::SiteType IOtype,
                                 geometry::LatticeData* iLatDat,
                                 int iBoundaryId);
          std::vector<int> GatherProcList(bool hasBoundary);

          void ReadParameters(geometry::LatticeData::SiteType IOtype);
          void allocate();

          void InitialiseBoundaryDensities();
          void UpdateBoundaryDensities();

          void InitialiseCosCycle(int i);
          void UpdateCosCycle(int i);

          void InitialiseFromFile(int i);

          void SortValuesFromFile(std::vector<double> &time, std::vector<double> &value);

          int nTotIOlets;
          // Number of IOlets and vector of their indices for communication purposes
          int nIOlets;
          std::vector<int> iolets;

          distribn_t *density_cycle;
          unsigned long *density_period;
          distribn_t *density;
          distribn_t *density_avg;
          distribn_t *density_amp;
          distribn_t *density_phs;
          distribn_t *density_min;
          distribn_t *density_max;

          std::string *filename;
          // Can just check equality with empty string, but despite its simplicity it still
          // confused me so I use a bool array initialised at construction
          bool *read_from_file;

          SimulationState* mState;
          SimConfig* mSimConfig;
          util::UnitConverter* mUnits;
      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_BOUNDARYVALUES_H */
