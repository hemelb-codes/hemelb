#ifndef HEMELB_LB_STREAMERS_LBGKNN_H
#define HEMELB_LB_STREAMERS_LBGKNN_H

#include "lb/kernels/BaseKernel.h"
#include "lb/SimulationState.h"
#include <cassert>
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      // Forward declaration needed by the struct
      template<class tRheologyModel> class LBGKNN;

      template<class tRheologyModel>
      struct HydroVars<LBGKNN<tRheologyModel> > : public HydroVarsBase
      {
        public:
          HydroVars(const distribn_t* const f) :
            HydroVarsBase(f)
          {

          }

          distribn_t tau;
      };

      /*
       * Class extending the original BGK collision operator to support non-Newtonian
       * fluids. Implements support for relaxation time not constant across the domain.
       */
      template<class tRheologyModel>
      class LBGKNN : public BaseKernel<LBGKNN<tRheologyModel> >
      {
        public:

          LBGKNN(InitParams& initParams)
          {
            InitState(initParams);
          }

          void DoCalculateDensityVelocityFeq(HydroVars<LBGKNN>& hydroVars, site_t index)
          {
            D3Q15::CalculateDensityVelocityFEq(hydroVars.f,
                                               hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z,
                                               hydroVars.f_eq);

            // Use the value of tau computed during the previous time step in coming calls to DoCollide
            assert(index < (site_t) mTau.size());
            hydroVars.tau = mTau[index];

            // Compute the local relaxation time that will be used in the next time step
            UpdateLocalTau(mTau[index], hydroVars);
          }

          void DoCalculateFeq(HydroVars<LBGKNN>& hydroVars, site_t index)
          {
            D3Q15::CalculateFeq(hydroVars.density,
                                hydroVars.v_x,
                                hydroVars.v_y,
                                hydroVars.v_z,
                                hydroVars.f_eq);

            // Use the value of tau computed during the previous time step in coming calls to DoCollide
            assert(index < (site_t) mTau.size());
            hydroVars.tau = mTau[index];

            // Compute the local relaxation time that will be used in the next time step
            UpdateLocalTau(mTau[index], hydroVars);
          }

          distribn_t DoCollide(const LbmParameters* const lbmParams,
                               HydroVars<LBGKNN>& hydroVars,
                               unsigned int direction)
          {
            double omega = -1.0 / hydroVars.tau;
            return hydroVars.f[direction] + hydroVars.f_neq[direction] * omega;
          }

          void DoReset(InitParams* initParams)
          {
            InitState(*initParams);
          }

          /*
           *  Helper method used in testing in order to access the mTau array after
           *  being set by DoCalculateDensityVelocityFeq
           */
          const std::vector<distribn_t>& GetTauValues() const
          {
            return mTau;
          }

        private:
          /*
           * Vector containing the current relaxation time for each site in the domain. It will be initialised
           * with the relaxation time corresponding to HemeLB's default Newtonian viscosity and each time step
           * will be updated based on the local hydrodynamic configuration
           */
          std::vector<distribn_t> mTau;

          /* Current time step */
          distribn_t mTimeStep;

          /* Current space step */
          distribn_t mSpaceStep;

          /*
           *  Helper method to set/update member variables. Called from the constructor and Reset()
           *
           *  @param initParams struct used to store variables required for initialisation of various operators
           */
          void InitState(const InitParams& initParams)
          {
            // Initialise relaxation time across the domain to HemeLB's default value.
            mTau.resize(initParams.siteCount, initParams.lbmParams->GetTau());
            mTimeStep = initParams.lbmParams->GetTimeStep();
            mSpaceStep = initParams.latDat->GetVoxelSize();
          }

          /*
           *  Helper method to update the value of local relaxation time (tau) from a given hydrodynamic
           *  configuration. It requires values of f_neq and density at the current time step and it will
           *  compute the value of tau to be used in the next time step.
           *
           *  @param localTau input: tau being used during the current time step,
           *                  output: tau to be used in the following time step
           *  @param hydroVars hydrodynamic configuration for a given lattice site
           */
          void UpdateLocalTau(distribn_t& localTau, HydroVars<LBGKNN>& hydroVars) const
          {

            /*
             *  TODO optimise, at this point hydroVars.f_neq has not been computed yet, so we
             *  need to do it here for the shear-rate calculator. However, the streamer will do
             *  it again before DoCollide is called.
             *
             *  Modify *all* the kernels to take care of this operation.
             */
            for (unsigned f_index = 0; f_index < D3Q15::NUMVECTORS; f_index++)
            {
              hydroVars.f_neq[f_index] = hydroVars.f[f_index] - hydroVars.f_eq[f_index];
            }

            /*
             * Shear-rate returned by CalculateShearRate is dimensionless and CalculateTauForShearRate
             * wants it in units of s^{-1}
             */
            double shear_rate = D3Q15::CalculateShearRate(localTau, hydroVars.f_neq, hydroVars.density) / mTimeStep;

            // Update tau
            localTau = tRheologyModel::CalculateTauForShearRate(shear_rate,
                                                                hydroVars.density,
                                                                mSpaceStep,
                                                                mTimeStep);

            // In some rheology models viscosity tends to infinity as shear rate goes to zero.
            assert( !std::isinf(localTau) );
            assert( !std::isnan(localTau) );
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LBGKNN_H */
