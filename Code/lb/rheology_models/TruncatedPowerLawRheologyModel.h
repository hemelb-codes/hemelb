#ifndef HEMELB_LB_RHEOLOGYMODELS_TRUNCATEDPOWERLAWRHEOLOGYMODEL_H_
#define HEMELB_LB_RHEOLOGYMODELS_TRUNCATEDPOWERLAWRHEOLOGYMODEL_H_

#include "lb/rheology_models/AbstractRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      static const double GAMMA_ZERO = 0.1937; // Pa^{1/2}
      static const double GAMMA_INF = 0.055; // (Pa s)^{1/2}

      static const double M_CONSTANT = 1e-3; // TODO units
      static const double N_CONSTANT = 0.5; // TODO units n<1 shear thinning, n>1 shear thickening

      class TruncatedPowerLawRheologyModel : public AbstractRheologyModel
      {
        public:
          static double CalculateTauForShearRate(const double &iShearRate, const distribn_t &iDensity);
      };
    }
  }
}

#endif /* HEMELB_LB_RHEOLOGYMODELS_TRUNCATEDPOWERLAWRHEOLOGYMODEL_H_ */
