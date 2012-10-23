// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/kernels/rheologyModels/CassonRheologyModel.h"
#include "util/utilityFunctions.h"
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        double CassonRheologyModel::CalculateViscosityForShearRate(const double &iShearRate,
                                                                   const distribn_t &iDensity)
        {
          double k0_k1_gamma = K0 + K1 * sqrt(iShearRate);
          double eta = (k0_k1_gamma * k0_k1_gamma) / iShearRate;

          // In the Casson rheology model, viscosity tends to infinity as shear rate goes to zero. Bound it.
          eta = hemelb::util::NumericalFunctions::min(eta, CASSON_MAX_VISCOSITY);

          // TODO Investigate whether we should be using BLOOD_DENSITY_Kg_per_m3*iDensity
          double nu = eta / BLOOD_DENSITY_Kg_per_m3;

          return nu;
        }
      }
    }
  }
}
