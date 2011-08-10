#include "lb/collisions/implementations/LBGKNN.h"
#include "lb/rheology_models/CassonRheologyModel.h"
#include "lb/rheology_models/TruncatedPowerLawRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {
        template<class tNonNewtonianModel>
        std::vector<double> LBGKNN<tNonNewtonianModel>::mTau;

        template<class tNonNewtonianModel>
        size_t LBGKNN<tNonNewtonianModel>::mCurrentTauIndex;

        template<class tNonNewtonianModel>
        distribn_t LBGKNN<tNonNewtonianModel>::mCurrentDensity;

        template<class tNonNewtonianModel>
        void LBGKNN<tNonNewtonianModel>::createTauArray(const site_t& iSize,
                                                        const double& iDefaultTau)
        {
            LBGKNN<tNonNewtonianModel>::mTau.resize(iSize,iDefaultTau);
        }

        template<class tNonNewtonianModel>
        void LBGKNN<tNonNewtonianModel>::getSiteValues(const distribn_t* f,
                                                       distribn_t &density,
                                                       distribn_t &v_x,
                                                       distribn_t &v_y,
                                                       distribn_t &v_z,
                                                       distribn_t* f_eq,
                                                       const site_t index)
        {
          mCurrentTauIndex = index;
          mCurrentDensity = density;

          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_eq);
        }

//        void LBGK::getBoundarySiteValues(const distribn_t* f,
//                                         const distribn_t &density,
//                                         const distribn_t &v_x,
//                                         const distribn_t &v_y,
//                                         const distribn_t &v_z,
//                                         distribn_t* f_eq,
//                                         const site_t index)
//        {
//          D3Q15::CalculateFeq(density, v_x, v_y, v_z, f_eq);
//        }
//
//        void LBGK::doPostCalculations(const distribn_t* f,
//                                      const geometry::LatticeData* bLatDat,
//                                      const site_t inde)
//        {
//
//        }

        template<class tNonNewtonianModel>
        distribn_t LBGKNN<tNonNewtonianModel>::getOperatorElement(distribn_t &f_i,
                                                                  distribn_t &f_neq_i,
                                                                  const LbmParameters* iLbmParams)
        {
            double old_tau_value = LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex];
            double shear_rate = D3Q15::CalculateShearRate(old_tau_value, f_neq_i);
            LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex] = tNonNewtonianModel::CalculateTauForShearRate(shear_rate, mCurrentDensity);;

            return (f_neq_i/LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex]);
        }


        // Explicit instantiation
        template class LBGKNN<hemelb::lb::rheology_models::CassonRheologyModel>;
        template class LBGKNN<hemelb::lb::rheology_models::TruncatedPowerLawRheologyModel>;
      }
    }
  }
}
