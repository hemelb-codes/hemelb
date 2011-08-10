#include "lb/rheology_models/CassonRheologyModel.h"

#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double CassonRheologyModel::CalculateTauForShearRate(const double &iShearRate,
                                                           const distribn_t &iDensity)
      {
        double k0_k1_gamma = K0 + K1*sqrt(iShearRate);
        double eta = (k0_k1_gamma*k0_k1_gamma)/iShearRate;

        return (3*eta)/iDensity + 1/2;
      }
    }
  }
}
