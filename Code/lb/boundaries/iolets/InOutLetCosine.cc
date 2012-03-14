#include "lb/boundaries/iolets/InOutLetCosine.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        InOutLetCosine::InOutLetCosine() :
          InOutLetCycle<1, false> ()
        {

        }

        void InOutLetCosine::DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig)
        {
          iSimConfig->DoIOForCosineInOutlet(iParent, iIsLoading, this);
        }

        InOutLet* InOutLetCosine::Clone()
        {
          InOutLetCosine* copy = new InOutLetCosine(*this);

          return copy;
        }

        InOutLetCosine::~InOutLetCosine()
        {

        }

        void InOutLetCosine::CalculateCycle(std::vector<distribn_t> &densityCycle, const SimulationState *iState)
        {
          double w = 2.0 * PI / (double) iState->GetTotalTimeSteps();

          for (unsigned int time_step = 0; time_step < densityCycle.size(); time_step++)
          {
            densityCycle[time_step] = GetDensityMean() + GetDensityAmp() * cos(w * (double) (time_step
                + iState->Get0IndexedTimeStep()) + Phase);
          }
        }

        distribn_t InOutLetCosine::GetDensityMean()
        {
          return  mUnits->ConvertPressureToLatticeUnits(PressureMeanPhysical) / Cs2;
        }

        distribn_t InOutLetCosine::GetDensityAmp()
        {
          return mUnits->ConvertPressureGradToLatticeUnits(PressureAmpPhysical) / Cs2;
        }

      }
    }
  }
}
