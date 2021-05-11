// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/kernels/rheologyModels/AbstractRheologyModel.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        template<class tRheologyImplementation>
        double AbstractRheologyModel<tRheologyImplementation>::CalculateTauForShearRate(
            const double &iShearRate, const distribn_t &iDensity, const LbmParameters& lbParams) const
        {
	  auto&& dx = lbParams.GetVoxelSize();
	  auto&& dt = lbParams.GetTimeStep();
	  auto self = static_cast<const tRheologyImplementation*>(this);
	  auto eta = self->CalculateViscosityForShearRate(iShearRate, iDensity);
	  // TODO: Investigate whether we should be using lbParams.GetFluidDensity() * iDensity
	  auto nu = eta / lbParams.GetFluidDensity();
          return 0.5 + (dt * nu) / (Cs2 * dx * dx);
        }

        // Explicit instantiation (a way of splitting templated classes into .h and .cc files)
        template class AbstractRheologyModel<CassonRheologyModel>;
        template class AbstractRheologyModel<TruncatedPowerLawRheologyModel>;
        template class AbstractRheologyModel<CarreauYasudaRheologyModel<HumanCYFit> >;
        template class AbstractRheologyModel<CarreauYasudaRheologyModel<MouseCYFit> >;
      }
    }
  }
}
