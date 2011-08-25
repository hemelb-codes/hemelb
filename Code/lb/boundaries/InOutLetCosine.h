#ifndef HEMELB_LB_BOUNDARIES_INOUTLETCOSINE_H
#define HEMELB_LB_BOUNDARIES_INOUTLETCOSINE_H

#include "lb/boundaries/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      class InOutLetCosine : public InOutLet
      {
        public:
          InOutLetCosine(unsigned long iPeriod,
                         double iPMin,
                         double iPMax,
                         util::Vector3D iPosition,
                         util::Vector3D iNormal,
                         util::UnitConverter* iUnits,
                         double iPMean,
                         double iPAmp,
                         double iPPhase);
          virtual ~InOutLetCosine();

          virtual void CalculateCycle(std::vector<distribn_t> &density_cycle,
                                      const SimulationState *iState);

          virtual void ResetValues();

        private:
          double PressureMeanPhysical;
          distribn_t DensityMeanLattice;

          double PressureAmpPhysical;
          distribn_t DensityAmpLattice;

          double Phase;

      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_INOUTLETCOSINE_H */
