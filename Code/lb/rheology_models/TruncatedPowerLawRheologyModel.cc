#include "lb/rheology_models/TruncatedPowerLawRheologyModel.h"

#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double TruncatedPowerLawRheologyModel::CalculateTauForShearRate(const double &iShearRate,
                                                           const distribn_t &iDensity)
      {

        double gamma = hemelb::util::NumericalFunctions::max(GAMMA_ZERO,
                                                             hemelb::util::NumericalFunctions::min(GAMMA_INF, iShearRate));
        double nu = M_CONSTANT * pow(gamma, N_CONSTANT-1);
        return (6*nu + 1)/2;
      }
    }
  }
}



