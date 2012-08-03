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

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        double CarreauYasudaRheologyModel::CalculateViscosityForShearRate(const double &iShearRate,
                                                                          const distribn_t &iDensity)
        {
          double eta = ETA_INF + (ETA_ZERO - ETA_INF) * pow( (1.0 + pow(LAMBDA * iShearRate, A)),
                                                             (N - 1.0) / A);

          // TODO Investigate whether we should be using BLOOD_DENSITY_Kg_per_m3*iDensity
          double nu = eta / BLOOD_DENSITY_Kg_per_m3;

          return nu;
        }
      }
    }
  }
}
