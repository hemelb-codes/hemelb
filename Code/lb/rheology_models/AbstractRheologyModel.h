#ifndef HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_
#define HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_

#include "constants.h"

namespace hemelb
{
  namespace lb
  {
    namespace rheology_models
    {
      class AbstractRheologyModel
      {
        public:
          /*
           *  Computes the local viscosity (eta) corresponding to shear rate, density,
           *  voxel size and time step duration according to a given rheology model.
           *
           *  We do not implement this method here and let subclasses do it according
           *  to different rheology models.
           *
           *  @param iShearRate local shear rate value.
           *  @param iDensity local density. TODO at the moment this value is not used
           *         in any subclass.
           *  @param iVoxelSize voxel size.
           *  @param iTimeStepsPerCycle number of time steps per cycle. Used to compute
           *         time step duration.
           */
          static double CalculateTauForShearRate(const double &iShearRate,
                                                 const distribn_t &iDensity,
                                                 const double &iVoxelSize,
                                                 const double &iTimeStep);

          /*
           *  Computes the relaxation parameter tau from a given kinematic viscosity nu
           *  according to equation (2.37) in Marco's thesis.
           *
           *  tau = 0.5 + (timestep * eta) / (Cs2 * density * voxelsize^2)
           *
           *  Cs2 is the dimensionless speed of the sound squared.
           *
           *  @param iNu kinematic viscosity (m^2/s)
           */
          static double CalculateTauForViscosity(const double &iNu,
                                                 const double &iTimeStep,
                                                 const double &iVoxelSize);
        private:
          AbstractRheologyModel();
      };
    }
  }
}


#endif /* HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_ */
