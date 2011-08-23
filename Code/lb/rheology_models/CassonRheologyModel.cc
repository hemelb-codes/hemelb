#include "lb/rheology_models/CassonRheologyModel.h"

#include <math.h>
#include <cassert>
#include <iostream>

#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      double CassonRheologyModel::CalculateTauForShearRate(const double &iShearRate,
                                                           const distribn_t &iDensity,
                                                           const double &iVoxelSize,
                                                           const unsigned &iTimeStepsPerCycle)
      {
        double k0_k1_gamma = K0 + K1*sqrt(iShearRate);
        double eta = (k0_k1_gamma*k0_k1_gamma)/iShearRate;

        // In the Casson rheology model, viscosity tends to infinity as shear rate goes to zero. Bound it.
        // TODO move this constant to the header file.
        eta = hemelb::util::NumericalFunctions::min(eta, 0.16);

        double timestep = PULSATILE_PERIOD_s / (double) iTimeStepsPerCycle;
        // TODO Investigate whether we should be using the default blood density or the value computed locally (iDensity)
        return 0.5 + (timestep * eta/BLOOD_DENSITY_Kg_per_m3) / (Cs2 * iVoxelSize * iVoxelSize);
      }
    }
  }
}
