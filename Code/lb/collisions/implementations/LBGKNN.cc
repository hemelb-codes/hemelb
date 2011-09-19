#include "lb/collisions/implementations/LBGKNN.h"
#include "lb/rheology_models/CassonRheologyModel.h"
#include "lb/rheology_models/TruncatedPowerLawRheologyModel.h"
#include "lb/rheology_models/CarreauYasudaRheologyModel.h"
#include <cmath>

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
        const geometry::LatticeData* LBGKNN<tNonNewtonianModel>::mLatDat;

        template<class tNonNewtonianModel>
        const SimulationState* LBGKNN<tNonNewtonianModel>::mState;


        template<class tNonNewtonianModel>
        void LBGKNN<tNonNewtonianModel>::createTauArray(const site_t& iSize,
                                                        const double& iDefaultTau)
        {
            LBGKNN<tNonNewtonianModel>::mTau.resize(iSize,iDefaultTau);
        }

        template<class tNonNewtonianModel>
        void LBGKNN<tNonNewtonianModel>::setStateObjects(const geometry::LatticeData* iLatDat,
                             const SimulationState* iState)
        {
            LBGKNN<tNonNewtonianModel>::mLatDat = iLatDat;
            LBGKNN<tNonNewtonianModel>::mState = iState;
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

          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_eq);

          // TODO redundant, f_neq is computed again when returning control to DoStreamAndCollide, it could be optimised.
          distribn_t f_neq[D3Q15::NUMVECTORS];
          for (unsigned f_index = 0; f_index < D3Q15::NUMVECTORS; f_index++)
          {
            f_neq[f_index] = f[f_index] - f_eq[f_index];
          }

          assert(mCurrentTauIndex < LBGKNN<tNonNewtonianModel>::mTau.size());
          double old_tau_value = LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex];

          // TODO no need to compute timestep duration for *every* site every time step
          double timestep = PULSATILE_PERIOD_s / (double) mState->GetTimeStepsPerCycle();
          double shear_rate = D3Q15::CalculateShearRate(old_tau_value, f_neq, timestep, mLatDat->GetVoxelSize(), density);
          LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex] = tNonNewtonianModel::CalculateTauForShearRate(shear_rate,
                                                                                                            density,
                                                                                                            mLatDat->GetVoxelSize(),
                                                                                                            timestep);
          // In some rheology models viscosity tends to infinity as shear rate goes to zero.
          assert( !std::isinf(LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex]) );
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
            double omega = -1.0/LBGKNN<tNonNewtonianModel>::mTau[mCurrentTauIndex];
            return (omega * f_neq_i);
        }


        // Explicit instantiation (a way of splitting templated classes into .h and .cc files)
        template class LBGKNN<hemelb::lb::rheology_models::CassonRheologyModel>;
        template class LBGKNN<hemelb::lb::rheology_models::TruncatedPowerLawRheologyModel>;
        template class LBGKNN<hemelb::lb::rheology_models::CarreauYasudaRheologyModel>;
      }
    }
  }
}
