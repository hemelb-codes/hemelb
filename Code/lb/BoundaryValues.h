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

        void InitialiseBoundaryDensities();

        int nTotInlets, nTotOutlets;

        std::vector<distribn_t> inlet_density_cycle, outlet_density_cycle;
        distribn_t *inlet_density_avg, *outlet_density_avg;
        distribn_t *inlet_density_amp, *outlet_density_amp;
        distribn_t *inlet_density_phs, *outlet_density_phs;

        SimulationState* mState;
        SimConfig* mSimConfig;
        util::UnitConverter* mUnits;
    };

  }
}

#endif /* HEMELB_LB_BOUNDARYVALUES_H */
