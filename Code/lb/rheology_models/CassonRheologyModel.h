#ifndef HEMELB_LB_RHEOLOGYMODELS_CASSONRHEOLOGYMODEL_H_
#define HEMELB_LB_RHEOLOGYMODELS_CASSONRHEOLOGYMODEL_H_

#include "lb/rheology_models/AbstractRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      static const double K0 = 0.1937; // Pa^{1/2}
      static const double K1 = 0.055; // (Pa s)^{1/2}

      class CassonRheologyModel : public AbstractRheologyModel
      {
        public:
          static double CalculateTauForShearRate(const double &iShearRate, const distribn_t &iDensity);
      };
    }
  }
}

#endif /* HEMELB_LB_RHEOLOGYMODELS_CASSONRHEOLOGYMODEL_H_ */
