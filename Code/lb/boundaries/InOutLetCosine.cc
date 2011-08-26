#include "lb/boundaries/InOutLetCosine.h"
#include "SimConfig.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      InOutLetCosine::InOutLetCosine() :
        InOutLetCycle<1, false>()
      {

      }

      void InOutLetCosine::DoIO(TiXmlElement *iParent, bool iIsLoading, SimConfig* iSimConfig)
      {
        iSimConfig->DoIO(iParent, iIsLoading, this);
      }

      InOutLet* InOutLetCosine::Clone()
      {
        InOutLetCosine* copy = new InOutLetCosine();
        copy->PressureMinPhysical = this->PressureMinPhysical;
        copy->PressureMaxPhysical = this->PressureMaxPhysical;
        copy->Position = this->Position;
        copy->Normal = this->Normal;
        copy->PressureMeanPhysical = this->PressureMeanPhysical;
        copy->PressureAmpPhysical = this->PressureAmpPhysical;
        copy->Phase = this->Phase;

        return copy;
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
              * (double) (time_step + iState->GetTimeStep() - 1) + Phase);
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

      distribn_t InOutLetCosine::GetDensityMean()
      {
        return DensityMeanLattice;
      }

      distribn_t InOutLetCosine::GetDensityAmp()
      {
        return DensityAmpLattice;
      }

    }
  }
}
