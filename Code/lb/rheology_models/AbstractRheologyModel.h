#ifndef HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_
#define HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_

#include "constants.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      class AbstractRheologyModel
      {
        public:
          static double CalculateTauForShearRate(const double &iShearRate,
                                                 const distribn_t &iDensity);

        private:
          AbstractRheologyModel();
      };
    }
  }
}


#endif /* HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_ */
