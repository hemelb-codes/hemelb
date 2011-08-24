#include "lb/rheology_models/AbstractRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double AbstractRheologyModel::CalculateTauForViscosity(const double &iNu,
                                                             const double &iTimeStep,
                                                             const double &iVoxelSize)
      {
        return 0.5 + (iTimeStep * iNu) / (Cs2 * iVoxelSize * iVoxelSize);
      }
    }
  }
}
