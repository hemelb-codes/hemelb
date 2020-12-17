// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/kernels/rheologyModels/CarreauYasudaRheologyModel.h"
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {

        template<class CYFIT>
        double CarreauYasudaRheologyModel<CYFIT>::CalculateViscosityForShearRate(
            const double &iShearRate, const distribn_t &iDensity) const
        {
	  constexpr double exponent = (CYFIT::N - 1.0) / CYFIT::A;
          double eta = CYFIT::ETA_INF
              + (CYFIT::ETA_ZERO - CYFIT::ETA_INF)
                  * pow( (1.0 + pow(CYFIT::LAMBDA * iShearRate, CYFIT::A)),
                        exponent);
	  return eta;
        }

        template class CarreauYasudaRheologyModel<HumanCYFit>;
        template class CarreauYasudaRheologyModel<MouseCYFit>;
      }
    }
  }
}
