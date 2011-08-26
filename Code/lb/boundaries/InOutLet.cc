#include "lb/boundaries/InOutLet.h"
#include "debug/Debugger.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      InOutLet::InOutLet()
      {

      }

      InOutLet::~InOutLet()
      {

      }

      void InOutLet::InitialiseCycle(std::vector<distribn_t> &density_cycle,
                                     const SimulationState *iState)
      {
        // Should be done anyway
        ResetValues();

        if (UpdatePeriod == 0)
        {
          density_cycle.resize(iState->GetTimeStepsPerCycle());
        }
        else
        {
          density_cycle.resize(UpdatePeriod);
        }

        CalculateCycle(density_cycle, iState);
      }

      void InOutLet::UpdateCycle(std::vector<distribn_t> &density_cycle,
                                 const SimulationState *iState)
      {
        if (UpdatePeriod == 0)
        {
          return;
        }
        else if ( (iState->GetTimeStep() - 1) % UpdatePeriod == 0)
        {
          CalculateCycle(density_cycle, iState);
        }
      }

      void InOutLet::ResetValues()
      {
        ResetCommonLatticeValues();
      }

      void InOutLet::ResetCommonLatticeValues()
      {
        DensityMinLattice = util::UnitConverter::ConvertPressureToLatticeUnits(PressureMinPhysical)
            / Cs2;
        DensityMaxLattice = util::UnitConverter::ConvertPressureToLatticeUnits(PressureMaxPhysical)
            / Cs2;
      }

      distribn_t InOutLet::GetDensityMin()
      {
        return DensityMinLattice;
      }

      distribn_t InOutLet::GetDensityMax()
      {
        return DensityMaxLattice;
      }

    }
  }
}
