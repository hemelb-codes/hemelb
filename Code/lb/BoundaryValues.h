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
        BoundaryValues(BoundaryComms* iInletComms,
                       BoundaryComms* iOutletComms,
                       geometry::LatticeData* iLatDat,
                       SimConfig* iSimConfig,
                       SimulationState* iSimState,
                       util::UnitConverter* iUnits);
        ~BoundaryValues();

        void CalculateBC(distribn_t f[],
                         hemelb::geometry::LatticeData::SiteType iSiteType,
                         unsigned int iBoundaryId,
                         distribn_t *density,
                         distribn_t *vx,
                         distribn_t *vy,
                         distribn_t *vz,
                         distribn_t f_neq[],
                         BoundaryComms* iInletComms,
                         BoundaryComms* iOutletComms);

        distribn_t GetInitialDensity();

        distribn_t GetInletDensityMin(int iBoundaryId);
        distribn_t GetInletDensityMax(int iBoundaryId);
        distribn_t GetOutletDensityMin(int iBoundaryId);
        distribn_t GetOutletDensityMax(int iBoundaryId);

        void Reset();

      private:
        void ReadParameters();
        void allocateInlets();
        void allocateOutlets();

        void FindIOletDensityExtrema();

        void InitialiseBoundaryDensities();
        void InitialiseCosCycle(int i,
                                int IOlets,
                                distribn_t* density_avg,
                                distribn_t* density_amp,
                                distribn_t* density_phs,
                                std::vector<distribn_t> &density_cycle);
        void InitialiseFromFile(int i,
                                std::string &filename,
                                std::vector<distribn_t> &density_cycle);

        int nTotInlets, nTotOutlets;

        std::vector<distribn_t> inlet_density_cycle, outlet_density_cycle;
        distribn_t *inlet_density_avg, *outlet_density_avg;
        distribn_t *inlet_density_amp, *outlet_density_amp;
        distribn_t *inlet_density_phs, *outlet_density_phs;
        distribn_t *inlet_density_min, *outlet_density_min;
        distribn_t *inlet_density_max, *outlet_density_max;

        std::string *inlet_file, *outlet_file;

        SimulationState* mState;
        SimConfig* mSimConfig;
        util::UnitConverter* mUnits;
    };

  }
}

#endif /* HEMELB_LB_BOUNDARYVALUES_H */
