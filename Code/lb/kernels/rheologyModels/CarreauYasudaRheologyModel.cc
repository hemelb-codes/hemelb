// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/kernels/rheologyModels/CarreauYasudaRheologyModel.h"
#include <cmath>

// Macro used to initialise struct NAME with a particular fit of the Carrea-Yasuda model
#define CY_FIT_INIT(NAME, ETA_INF_PA_S, ETA_ZERO_PA_S, LAMBDA_S, A_DIMLESS, N_DIMLESS) \
		    const double NAME::ETA_INF = ETA_INF_PA_S; \
        const double NAME::ETA_ZERO = ETA_ZERO_PA_S; \
        const double NAME::LAMBDA = LAMBDA_S; \
        const double NAME::A = A_DIMLESS; \
        const double NAME::N = N_DIMLESS

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        CY_FIT_INIT(HumanCYFit, 0.0035, 0.16, 8.2, 0.64, 0.2128);
        CY_FIT_INIT(MouseCYFit, 3.265e-3, 14.49e-3, 0.1829, 2.707, 0.4136);

        template<class CYFIT>
        double CarreauYasudaRheologyModel<CYFIT>::CalculateViscosityForShearRate(
            const double &iShearRate, const distribn_t &iDensity)
        {
          double eta = CYFIT::ETA_INF
              + (CYFIT::ETA_ZERO - CYFIT::ETA_INF)
                  * pow( (1.0 + pow(CYFIT::LAMBDA * iShearRate, CYFIT::A)),
                        (CYFIT::N - 1.0) / CYFIT::A);

          // TODO Investigate whether we should be using BLOOD_DENSITY_Kg_per_m3*iDensity
          double nu = eta / BLOOD_DENSITY_Kg_per_m3;

          return nu;
        }

        template class CarreauYasudaRheologyModel<HumanCYFit> ;
        template class CarreauYasudaRheologyModel<MouseCYFit> ;
      }
    }
  }
}
