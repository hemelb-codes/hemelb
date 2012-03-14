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
        :comms(NULL)
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
          return mUnits->ConvertPressureToLatticeUnits(GetPressureMin()) / Cs2;
        }

        distribn_t InOutLet::GetDensityMax()
        {
          return mUnits->ConvertPressureToLatticeUnits(GetPressureMax()) / Cs2;
        }

      }
    }
  }
}
