#include "lb/kernels/rheologyModels/AbstractRheologyModel.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace rheologyModels
      {
        template<class tRheologyImplementation>
        double AbstractRheologyModel<tRheologyImplementation>::CalculateTauForShearRate(const double &iShearRate,
                                                                               const distribn_t &iDensity,
                                                                               const double &iVoxelSize,
                                                                               const double &iTimeStep)
        {
          double nu = tRheologyImplementation::CalculateViscosityForShearRate(iShearRate, iDensity);
          return 0.5 + (iTimeStep * nu) / (Cs2 * iVoxelSize * iVoxelSize);
        }

        // Explicit instantiation (a way of splitting templated classes into .h and .cc files)
        template class AbstractRheologyModel<CassonRheologyModel> ;
        template class AbstractRheologyModel<TruncatedPowerLawRheologyModel> ;
        template class AbstractRheologyModel<CarreauYasudaRheologyModel> ;
      }
    }
  }
}
