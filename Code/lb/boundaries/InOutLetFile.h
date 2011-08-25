#ifndef HEMELB_LB_BOUNDARIES_INOUTLETFILE_H
#define HEMELB_LB_BOUNDARIES_INOUTLETFILE_H

#include "lb/boundaries/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      class InOutLetFile : public InOutLet
      {
        public:
          InOutLetFile(unsigned long iPeriod,
                       double iPMin,
                       double iPMax,
                       util::Vector3D iPosition,
                       util::Vector3D iNormal,
                       std::string iPFilePath);
          virtual ~InOutLetFile();

          virtual void CalculateCycle(std::vector<distribn_t> &density_cycle,
                                   const SimulationState *iState);

        private:
          std::string PressureFilePath;
      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_INOUTLETFILE_H */
