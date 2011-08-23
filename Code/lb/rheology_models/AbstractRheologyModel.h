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
           *  The viscosity value is then converted into the corresponding relaxation
           *  value tau according to equation ( TODO ) in Marco's thesis.
           *
           *  tau = 0.5 + (timestep * eta) / (Cs2 * density * voxelsize^2)
           *
           *  Cs2 is the dimensionless speed of the sound squared.
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
                                                 const unsigned &iTimeStepsPerCycle);
        private:
          AbstractRheologyModel();
      };
    }
  }
}


#endif /* HEMELB_LB_RHEOLOGYMODELS_ABSTRACTRHEOLOGYMODEL_H_ */
