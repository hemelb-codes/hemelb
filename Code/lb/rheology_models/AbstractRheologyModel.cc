#include "lb/rheology_models/AbstractRheologyModel.h"
#include "lb/rheology_models/RheologyModels.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      template <class CONCRETE_CLASS>
      double AbstractRheologyModel<CONCRETE_CLASS>::CalculateTauForShearRate(const double &iShearRate,
                                                                             const distribn_t &iDensity,
                                                                             const double &iVoxelSize,
                                                                             const double &iTimeStep)
      {
        double nu = CONCRETE_CLASS::CalculateViscosityForShearRate(iShearRate, iDensity);
        return 0.5 + (iTimeStep * nu) / (Cs2 * iVoxelSize * iVoxelSize);
      }

      // Explicit instantiation (a way of splitting templated classes into .h and .cc files)
      template class AbstractRheologyModel<CassonRheologyModel>;
      template class AbstractRheologyModel<TruncatedPowerLawRheologyModel>;
      template class AbstractRheologyModel<CarreauYasudaRheologyModel>;
    }
  }
}
