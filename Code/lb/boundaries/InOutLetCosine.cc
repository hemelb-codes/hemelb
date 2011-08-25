#include "lb/boundaries/InOutLetCosine.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      InOutLetCosine::InOutLetCosine(unsigned long iPeriod,
                                     double iPMin,
                                     double iPMax,
                                     util::Vector3D iPosition,
                                     util::Vector3D iNormal,
                                     double iPMean,
                                     double iPAmp,
                                     double iPPhase) :
        InOutLet(iPeriod, iPMin, iPMax, iPosition, iNormal), PressureMeanPhysical(iPMean),
            PressureAmpPhysical(iPAmp), Phase(iPPhase)
      {

      }

      InOutLetCosine::~InOutLetCosine()
      {

      }

      void InOutLetCosine::CalculateCycle(std::vector<distribn_t> &density_cycle,
                                          const SimulationState *iState)
      {
        double w = 2.0 * PI / (double) iState->GetTimeStepsPerCycle();

        for (unsigned int time_step = 0; time_step < density_cycle.size(); time_step++)
        {
          density_cycle[time_step] = DensityMeanLattice + DensityAmpLattice * cos(w
              * (double) time_step + Phase);
        }
      }

      void InOutLetCosine::ResetValues()
      {
        DensityMeanLattice
            = util::UnitConverter::ConvertPressureToLatticeUnits(PressureMeanPhysical) / Cs2;
        DensityAmpLattice
            = util::UnitConverter::ConvertPressureGradToLatticeUnits(PressureAmpPhysical) / Cs2;

        ResetCommonLatticeValues();
      }
    }
  }
}
