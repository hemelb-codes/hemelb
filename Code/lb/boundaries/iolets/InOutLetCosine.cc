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
          InOutLet(),PressureMeanPhysical(0.0),PressureAmpPhysical(0.0),Phase(0.0),Period(1.0)
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

        LatticeDensity InOutLetCosine::GetDensity(unsigned long time_step)
        {
          double w = 2.0 * PI / Period;
          return GetDensityMean() + GetDensityAmp() * cos(w * mUnits->ConvertTimeStepToPhysicalUnits(time_step) + Phase);
        }

        LatticeDensity InOutLetCosine::GetDensityMean()
        {
          return  mUnits->ConvertPressureToLatticeUnits(PressureMeanPhysical) / Cs2;
        }

        LatticeDensity InOutLetCosine::GetDensityAmp()
        {
          return mUnits->ConvertPressureGradToLatticeUnits(PressureAmpPhysical) / Cs2;
        }

      }
    }
  }
}
