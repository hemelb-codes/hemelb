// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H
#define HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H

#include "lb/kernels/rheologyModels/AbstractRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        // Carreau-Yasuda model constants
        static const double ETA_INF = 0.0035; // Pa*s
        static const double ETA_ZERO = 0.16; // Pa*s
        static const double LAMBDA = 8.2; // s
        static const double A = 0.64; //  dimensionless
        static const double N = 0.2128; //  dimensionless

        class CarreauYasudaRheologyModel : public AbstractRheologyModel<CarreauYasudaRheologyModel>
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
      }
    }
  }
}

#endif /* HEMELB_LB_KERNELS_RHEOLOGYMODELS_CARREAUYASUDARHEOLOGYMODEL_H */
