#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCOSINE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCOSINE_H

#include "lb/boundaries/iolets/InOutLetCycle.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        class InOutLetCosine : public InOutLetCycle<1, false>
        {
          public:
            InOutLetCosine();
            virtual ~InOutLetCosine();

            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, SimConfig* iSimConfig);
            virtual InOutLet* Clone();

            virtual void CalculateCycle(std::vector<distribn_t> &density_cycle,
                                        const SimulationState *iState);

            virtual void ResetValues();

            distribn_t GetDensityMean();
            distribn_t GetDensityAmp();

            double PressureMeanPhysical;
            double PressureAmpPhysical;

            double Phase;

          private:
            distribn_t DensityMeanLattice;
            distribn_t DensityAmpLattice;

        };

      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETCOSINE_H */
