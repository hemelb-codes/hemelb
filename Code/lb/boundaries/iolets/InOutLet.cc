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

        void InOutLet::Initialise(const util::UnitConverter* units)
        {
          mUnits = units;
        }

        distribn_t InOutLet::GetDensityMin()
        {
          return mUnits->ConvertPressureToLatticeUnits(PressureMinPhysical) / Cs2;
        }

        distribn_t InOutLet::GetDensityMax()
        {
          return mUnits->ConvertPressureToLatticeUnits(PressureMaxPhysical) / Cs2;
        }

      }
    }
  }
}
