#ifndef HEMELB_LB_BOUNDARIES_INOUTLET_H
#define HEMELB_LB_BOUNDARIES_INOUTLET_H

#include "util/Vector3D.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      class InOutLet
      {
        public:
          InOutLet(unsigned long iPeriod,
                   double iPMin,
                   double iPMax,
                   util::Vector3D iPosition,
                   util::Vector3D iNormal,
                   util::UnitConverter* iUnits);
          virtual ~InOutLet();

          // Should be called before simulation starts running (including after a reset)
          // Resizes density_cycle and calls CalculateCycle
          void InitialiseCycle(std::vector<distribn_t> &density_cycle,
                               const SimulationState *iState);
          void UpdateCycle(std::vector<distribn_t> &density_cycle, const SimulationState *iState);

          virtual void CalculateCycle(std::vector<distribn_t> &density_cycle,
                                      const SimulationState *iState) = 0;

          virtual void ResetValues();
          void ResetCommonLatticeValues();

          distribn_t density;

          distribn_t GetDensityMin();
          distribn_t GetDensityMax();

          util::Vector3D GetPosition();
          util::Vector3D GetNormal();

        protected:
          util::UnitConverter* mUnits;
          const unsigned long UpdatePeriod;

        private:
          const double PressureMinPhysical;
          distribn_t DensityMinLattice;

          const double PressureMaxPhysical;
          distribn_t DensityMaxLattice;

          const util::Vector3D Position;
          const util::Vector3D Normal;

      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_INOUTLET_H */
