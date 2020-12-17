// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H
#define HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H

#include "lb/kernels/rheologyModels/AbstractRheologyModel.h"

// Macro used to define a fit of the Carreau-Yasuda model. We chose this design
// (with multiple structs containing the same static const variables instead of
// inheritance or similar) because CalculateViscosityForShearRate is performance
// critical and we want as much arithmetic done at compile time as possible.

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {

        template<class CYFIT>
        class CarreauYasudaRheologyModel : public AbstractRheologyModel<
            CarreauYasudaRheologyModel<CYFIT> >
        {
          public:
            CarreauYasudaRheologyModel(const InitParams&) {}
            /*
             *  Compute nu for a given shear rate according to the Carreau-Yasuda model:
             *
             *  eta = ETA_INF + (ETA_ZERO - ETA_INF) * (1 + (LAMBDA*iShearRate)^A)^((N-1)/A)
             *  nu = eta / density
             *
             *  @param iShearRate local shear rate value (s^{-1}).
             *  @param iDensity local density. TODO at the moment this value is not used
             *         in any subclass.
             *
             *  @return dynamic viscosity (Pa s).
             */
            double CalculateViscosityForShearRate(const double &iShearRate,
						  const distribn_t &iDensity) const;
        };

#define CY_FIT_NEW(NAME, ETA_INF_PA_S, ETA_ZERO_PA_S, LAMBDA_S, A_DIMLESS, N_DIMLESS) \
	struct NAME							\
	{								\
	  static constexpr double ETA_INF = ETA_INF_PA_S;		\
	  static constexpr double ETA_ZERO = ETA_ZERO_PA_S;		\
	  static constexpr double LAMBDA = LAMBDA_S;			\
	  static constexpr double A = A_DIMLESS;			\
	  static constexpr double N = N_DIMLESS;			\
	}

        CY_FIT_NEW(HumanCYFit, 0.0035, 0.16, 8.2, 0.64, 0.2128);
        CY_FIT_NEW(MouseCYFit, 3.265e-3, 14.49e-3, 0.1829, 2.707, 0.4136);
#undef CY_FIT_NEW

        using CarreauYasudaRheologyModelHumanFit = CarreauYasudaRheologyModel<HumanCYFit>;
        using CarreauYasudaRheologyModelMouseFit = CarreauYasudaRheologyModel<MouseCYFit>;
      }
    }
  }
}

#endif /* HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H */
