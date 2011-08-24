#ifndef HEMELB_LB_RHEOLOGYMODELS_CASSONRHEOLOGYMODEL_H_
#define HEMELB_LB_RHEOLOGYMODELS_CASSONRHEOLOGYMODEL_H_

#include "lb/rheology_models/AbstractRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      // Casson model constants
      static const double K0 = 0.1937; // Pa^{1/2}
      static const double K1 = 0.055; // (Pa*s)^{1/2}

      static const double CASSON_MAX_VISCOSITY = 0.16; // Pa*s

      class CassonRheologyModel : public AbstractRheologyModel
      {
        public:
          /*
           * Compute tau for a given shear rate according to the Casson model:
           *
           * eta = (K0 + K1*sqrt(iShearRate))^2 / iShearRate
           *
           * See AbstractRheologyModel.h for description of the arguments and how tau
           * is computed from eta.
           */
          static double CalculateTauForShearRate(const double &iShearRate,
                                                 const distribn_t &iDensity,
                                                 const double &iVoxelSize,
                                                 const double &iTimeStep);
      };
    }
  }
}

#endif /* HEMELB_LB_RHEOLOGYMODELS_CASSONRHEOLOGYMODEL_H_ */
