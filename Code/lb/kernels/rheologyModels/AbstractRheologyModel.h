
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H
#define HEMELB_LB_KERNELS_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H

#include "constants.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        template<class tRheologyImplementation>
        class AbstractRheologyModel
        {
          public:
            /*
             *  Computes the relaxation parameter tau according to equation (2.37) in
             *  Marco's thesis.
             *
             *  tau = 0.5 + (timestep * nu) / (Cs2 * voxelsize^2)
             *
             *  Cs2 is the dimensionless speed of the sound squared.
             *  nu is the kinematic viscosity (m^2/s).
             *
             *  This method relies on CalculateViscosityForShearRate to compute nu based
             *  on a given shear rate, density, and a given rheology model.
             *
             *  @param iShearRate local shear rate value (s^{-1}).
             *  @param iDensity local density. TODO at the moment this value is not used
             *         in any subclass.
             *  @param iVoxelSize voxel size.
             *  @param iTimeStep time step duration.
             *
             *  @return relaxation time (dimensionless)
             */
            static double CalculateTauForShearRate(const double &iShearRate,
                                                   const distribn_t &iDensity,
                                                   const double &iVoxelSize,
                                                   const double &iTimeStep);

            /*
             *  Computes the kinematic viscosity predicted by a given rheology model for
             *  given shear rate and viscosity
             *
             *  To be implemented in each AbstractRheologyModel subclass.
             *
             *  @param iShearRate local shear rate value (s^{-1}).
             *  @param iDensity local density. TODO at the moment this value is not used
             *         in any subclass.
             *
             *  @return kinematic viscosity (m^2/s).
             */
            static double CalculateViscosityForShearRate(const double &iShearRate,
                                                         const distribn_t &iDensity);
          private:
            AbstractRheologyModel();
        };
      }
    }
  }
}

#endif /* HEMELB_LB_KERNELS_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H */
