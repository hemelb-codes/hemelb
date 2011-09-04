#include "lb/rheology_models/TruncatedPowerLawRheologyModel.h"

#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double TruncatedPowerLawRheologyModel::CalculateViscosityForShearRate(const double &iShearRate,
                                                                            const distribn_t &iDensity)
      {
        // Don't allow shear rates outside [GAMMA_ZERO, GAMMA_INF]
        double gamma = hemelb::util::NumericalFunctions::max(GAMMA_ZERO,
                                                             hemelb::util::NumericalFunctions::min(GAMMA_INF, iShearRate));
        double nu = M_CONSTANT * pow(gamma, N_CONSTANT-1);
        return nu;
      }
    }
  }
}



