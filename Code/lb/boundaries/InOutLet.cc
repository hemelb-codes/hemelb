#include "lb/boundaries/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      InOutLet::InOutLet(unsigned long iPeriod,
                         double iPMin,
                         double iPMax,
                         util::Vector3D iPosition,
                         util::Vector3D iNormal) :
        UpdatePeriod(iPeriod), PressureMinPhysical(iPMin), PressureMaxPhysical(iPMax),
            Position(iPosition), Normal(iNormal)
      {

      }

      InOutLet::~InOutLet()
      {

      }

      void InOutLet::InitialiseCycle(std::vector<distribn_t> &density_cycle,
                                     const SimulationState *iState)
      {
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
        else
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

      util::Vector3D InOutLet::GetPosition()
      {
        return Position;
      }

      util::Vector3D InOutLet::GetNormal()
      {
        return Normal;
      }

    }
  }
}
