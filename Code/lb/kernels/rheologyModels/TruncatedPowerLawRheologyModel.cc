// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/kernels/BaseKernel.h"
#include "lb/kernels/rheologyModels/TruncatedPowerLawRheologyModel.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
	TruncatedPowerLawRheologyModel::TruncatedPowerLawRheologyModel(InitParams& initParams) : M_CONSTANT(initParams.lbmParams->GetEta()) {
	}
        double TruncatedPowerLawRheologyModel::CalculateViscosityForShearRate(
            const double &iShearRate, const distribn_t &iDensity) const
        {
          // Don't allow shear rates outside [GAMMA_ZERO, GAMMA_INF]
          double gamma = util::NumericalFunctions::enforceBounds(iShearRate, GAMMA_ZERO, GAMMA_INF);
          return M_CONSTANT * pow(gamma, N_CONSTANT - 1);
        }
      }
    }
  }
}

