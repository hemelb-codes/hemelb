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
                                     util::UnitConverter* iUnits,
                                     double iPMean,
                                     double iPAmp,
                                     double iPPhase) :
        InOutLet(iPeriod, iPMin, iPMax, iPosition, iNormal, iUnits), PressureMeanPhysical(iPMean),
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
        DensityMeanLattice = mUnits->ConvertPressureToLatticeUnits(PressureMeanPhysical) / Cs2;
        DensityAmpLattice = mUnits->ConvertPressureGradToLatticeUnits(PressureAmpPhysical) / Cs2;

        ResetCommonLatticeValues();
      }
    }
  }
}
