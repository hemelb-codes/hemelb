#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H

#include "util/Vector3D.h"
#include "util/UnitConverter.h"
#include "xml/tinyxml.h"

namespace hemelb
{

  class SimConfig;

  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        class InOutLet
        {
          public:
            InOutLet();
            virtual ~InOutLet();

            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, SimConfig* iSimConfig) = 0;

            virtual InOutLet* Clone() = 0;

            // Should be called before simulation starts running (including after a reset)
            // Resizes density_cycle and calls CalculateCycle
            virtual void InitialiseCycle(std::vector<distribn_t> &density_cycle,
                                         const SimulationState *iState) = 0;
            virtual void UpdateCycle(std::vector<distribn_t> &density_cycle,
                                     const SimulationState *iState) = 0;
            virtual void CalculateCycle(std::vector<distribn_t> &density_cycle,
                                        const SimulationState *iState) = 0;
            virtual bool DoComms() = 0;

            virtual void ResetValues();
            void ResetCommonLatticeValues();

            distribn_t GetDensityMin();
            distribn_t GetDensityMax();

            distribn_t density;
            unsigned long UpdatePeriod;

            double PressureMinPhysical;
            double PressureMaxPhysical;

            util::Vector3D Position;
            util::Vector3D Normal;

          private:
            distribn_t DensityMinLattice;
            distribn_t DensityMaxLattice;

        };

      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H */
