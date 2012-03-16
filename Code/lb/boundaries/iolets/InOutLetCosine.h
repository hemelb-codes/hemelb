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
            virtual InOutLet* Clone() const;
            virtual void Reset(SimulationState &state){
              //pass;
            }
            LatticeDensity GetDensityMean() const;
            LatticeDensity GetDensityAmp() const;
            LatticeDensity GetDensity(unsigned long time_step) const;
            PhysicalPressure GetPressureMin() const
            {
              return PressureMeanPhysical - PressureAmpPhysical;
            }
            PhysicalPressure GetPressureMax() const
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
