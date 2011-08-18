#ifndef HEMELB_LB_BOUNDARYVALUES_H
#define HEMELB_LB_BOUNDARYVALUES_H

#include "lb/BoundaryComms.h"

namespace hemelb
{
  namespace lb
  {

    class BoundaryValues
    {
      public:
        BoundaryValues(BoundaryComms* iComms,
                       geometry::LatticeData::SiteType IOtype,
                       geometry::LatticeData* iLatDat,
                       SimConfig* iSimConfig,
                       SimulationState* iSimState,
                       util::UnitConverter* iUnits);
        ~BoundaryValues();

        distribn_t GetDensityMin(int iBoundaryId);
        distribn_t GetDensityMax(int iBoundaryId);

        void ResetPrePeriodDoubling();
        void ResetPostPeriodDoubling();

        bool IsCurrentProcTheBCProc();

      private:
        proc_t BCproc;

        void ReadParameters(geometry::LatticeData::SiteType IOtype);
        void allocate();

        void FindDensityExtrema();

        void InitialiseBoundaryDensities();
        void InitialiseCosCycle(int i);
        void InitialiseFromFile(int i);
        void SortValuesFromFile(std::vector<double> &time, std::vector<double> &value);

        int nTotIOlets;

        std::vector<distribn_t> density_cycle;
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

#endif /* HEMELB_LB_BOUNDARYVALUES_H */
