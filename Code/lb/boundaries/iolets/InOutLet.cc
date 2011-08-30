#include "lb/boundaries/iolets/InOutLet.h"
#include "debug/Debugger.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        InOutLet::InOutLet()
        {

        }

        InOutLet::~InOutLet()
        {

        }

        void InOutLet::ResetValues()
        {
          ResetCommonLatticeValues();
        }

        void InOutLet::ResetCommonLatticeValues()
        {
          DensityMinLattice
              = util::UnitConverter::ConvertPressureToLatticeUnits(PressureMinPhysical) / Cs2;
          DensityMaxLattice
              = util::UnitConverter::ConvertPressureToLatticeUnits(PressureMaxPhysical) / Cs2;
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
}
