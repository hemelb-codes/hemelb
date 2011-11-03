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
            InOutLetCycle<1, false>()
        {

        }

        void InOutLetCosine::DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig)
        {
          iSimConfig->DoIO(iParent, iIsLoading, this);
        }

        InOutLet* InOutLetCosine::Clone()
        {
          InOutLetCosine* copy = new InOutLetCosine(*this);

          return copy;
        }

        InOutLetCosine::~InOutLetCosine()
        {

        }

        void InOutLetCosine::CalculateCycle(std::vector<distribn_t> &densityCycle,
                                            const SimulationState *iState)
        {
          double w = 2.0 * PI / (double) iState->GetTimeStepsPerCycle();

          for (unsigned int time_step = 0; time_step < densityCycle.size(); time_step++)
          {
            densityCycle[time_step] = DensityMeanLattice
                + DensityAmpLattice
                    * cos(w * (double) (time_step + iState->Get0IndexedTimeStep()) + Phase);
          }
        }

        void InOutLetCosine::ResetValues()
        {
          DensityMeanLattice = mUnits->ConvertPressureToLatticeUnits(PressureMeanPhysical) / Cs2;
          DensityAmpLattice = mUnits->ConvertPressureGradToLatticeUnits(PressureAmpPhysical) / Cs2;

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
}
