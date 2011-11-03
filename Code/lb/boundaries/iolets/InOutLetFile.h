#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETFILE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETFILE_H

#include "lb/boundaries/iolets/InOutLetCycle.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        /*
         * Template values are chosen to be tUpdatePeriod = 0 and tComms = false
         * If a trace is read from a file it should be done once and then stored
         * on each relevant proc. If memory is a concern tComms can be set to true
         * and then only the BCproc will keep the entire trace in memory
         * WARNING: - be cautious of setting tUpdatePeriod to something else other than
         * zero, because it may not be what you expect - see comments on CalculateCycle in
         * cc file.
         */
        class InOutLetFile : public InOutLetCycle<0, false>
        {
          public:
            InOutLetFile();
            virtual ~InOutLetFile();

            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig);
            virtual InOutLet* Clone();

            virtual void CalculateCycle(std::vector<distribn_t> &densityCycle,
                                        const SimulationState *iState);

            std::string PressureFilePath;
        };

      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETFILE_H */
