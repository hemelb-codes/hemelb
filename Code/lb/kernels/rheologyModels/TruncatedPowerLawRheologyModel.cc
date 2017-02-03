// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
        double TruncatedPowerLawRheologyModel::CalculateViscosityForShearRate(
            const double &iShearRate, const distribn_t &iDensity)
        {
          // Don't allow shear rates outside [GAMMA_ZERO, GAMMA_INF]
          double gamma = util::NumericalFunctions::enforceBounds(iShearRate, GAMMA_ZERO, GAMMA_INF);
          double eta = M_CONSTANT * pow(gamma, N_CONSTANT - 1);

          // TODO Investigate whether we should be using BLOOD_DENSITY_Kg_per_m3*iDensity
          double nu = eta / BLOOD_DENSITY_Kg_per_m3;

          return nu;
        }
      }
    }
  }
}

