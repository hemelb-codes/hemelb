#include "lb/rheology_models/CassonRheologyModel.h"
#include "util/utilityFunctions.h"
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double CassonRheologyModel::CalculateViscosityForShearRate(const double &iShearRate,
                                                                 const distribn_t &iDensity)
      {
        double k0_k1_gamma = K0 + K1*sqrt(iShearRate);
        double eta = (k0_k1_gamma*k0_k1_gamma)/iShearRate;

        // In the Casson rheology model, viscosity tends to infinity as shear rate goes to zero. Bound it.
        eta = hemelb::util::NumericalFunctions::min(eta, CASSON_MAX_VISCOSITY);

        // TODO Investigate whether we should be using BLOOD_DENSITY_Kg_per_m3*iDensity
        double nu = eta/BLOOD_DENSITY_Kg_per_m3;

        return nu;
      }
    }
  }
}
