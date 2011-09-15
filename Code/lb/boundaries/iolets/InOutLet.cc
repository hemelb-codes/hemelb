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
              = mUnits->ConvertPressureToLatticeUnits(PressureMinPhysical) / Cs2;
          DensityMaxLattice
              = mUnits->ConvertPressureToLatticeUnits(PressureMaxPhysical) / Cs2;
        }

        void InOutLet::Initialise(const util::UnitConverter* units)
        {
          mUnits = units;
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
