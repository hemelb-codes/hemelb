#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCOSINE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCOSINE_H

#include "lb/boundaries/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        /*
         * Template values are chosen to be tUpdatePeriod = 1 and tComms = false
         * The cosine pressure trace can be easily calculated locally by any proc and
         * there is no need to store values for any time step beyond the current one
         */
        class InOutLetCosine : public InOutLet
        {
          public:
            InOutLetCosine();
            virtual ~InOutLetCosine();
            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig);
            virtual InOutLet* Clone();
            virtual void Reset(SimulationState &state)
            {
              Period=state.GetTotalTimeSteps();
            }
            LatticeDensity GetDensityMean();
            LatticeDensity GetDensityAmp();
            LatticeDensity GetDensity(unsigned long time_step);
            PhysicalPressure GetPressureMin()
            {
              return PressureMeanPhysical - PressureAmpPhysical;
            }
            PhysicalPressure GetPressureMax()
            {
              return PressureMeanPhysical + PressureAmpPhysical;
            }

            PhysicalPressure PressureMeanPhysical;
            PhysicalPressure PressureAmpPhysical;

            double Phase;
            double Period;


        };

      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCOSINE_H */
