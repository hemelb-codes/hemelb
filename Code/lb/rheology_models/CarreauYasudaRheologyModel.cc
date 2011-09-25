#include "lb/rheology_models/CarreauYasudaRheologyModel.h"
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double CarreauYasudaRheologyModel::CalculateViscosityForShearRate(const double &iShearRate,
                                                                        const distribn_t &iDensity)
      {
        double  eta = ETA_INF + (ETA_ZERO - ETA_INF) * pow((1.0 + pow(LAMBDA*iShearRate,A)), (N-1.0)/A);

        // TODO Investigate whether we should be using BLOOD_DENSITY_Kg_per_m3*iDensity
        double nu = eta/BLOOD_DENSITY_Kg_per_m3;

        return nu;
      }
    }
  }
}
