
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
#define CY_FIT_NEW(NAME) \
		    struct NAME \
        { \
            static const double ETA_INF; \
            static const double ETA_ZERO; \
            static const double LAMBDA; \
            static const double A; \
            static const double N; \
        } \

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        CY_FIT_NEW(HumanCYFit);
        CY_FIT_NEW(MouseCYFit);

        template<class CYFIT>
        class CarreauYasudaRheologyModel : public AbstractRheologyModel<
            CarreauYasudaRheologyModel<CYFIT> >
        {
          public:
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
             *  @return kinematic viscosity (m^2/s).
             */
            static double CalculateViscosityForShearRate(const double &iShearRate,
                                                         const distribn_t &iDensity);
        };

        typedef CarreauYasudaRheologyModel<HumanCYFit> CarreauYasudaRheologyModelHumanFit;
        typedef CarreauYasudaRheologyModel<MouseCYFit> CarreauYasudaRheologyModelMouseFit;
      }
    }
  }
}

#endif /* HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H */
