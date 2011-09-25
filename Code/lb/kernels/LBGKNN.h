#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGKNN_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGKNN_H

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
      template<class tNonNewtonianModel> class LBGKNN;

      template<class tNonNewtonianModel>
      struct HydroVars<LBGKNN<tNonNewtonianModel> > : public HydroVarsBase
      {
        public:
          HydroVars(const distribn_t* const f) :
            HydroVarsBase(f)
          {

          }

          size_t index;
      };



      /*
       * Class extending the original BGK collision operator to support non-newtonian
       * fluids. Implements support for relaxation time not constant across the domain.
       */
      template<class tNonNewtonianModel>
      class LBGKNN : public BaseKernel< LBGKNN<tNonNewtonianModel> >
      {
        private:
          /*
           *  Helper method to set/update member variables. Called from the constructor and Reset()
           *
           *  @param initParams struct used to store variables required for initialisation of various operators
           */
          void InitState(InitParams& initParams)
          {
            mTau.resize(initParams.siteCount,initParams.lbmParams->Tau());
            mTimeStep = initParams.lbmParams->GetTimeStep();
            mSpaceStep = initParams.latDat->GetVoxelSize();
          }

        public:

          LBGKNN(InitParams& initParams)
          {
            InitState(initParams);
          }

          void DoCalculateDensityVelocityFeq(HydroVars<LBGKNN>& hydroVars, site_t index)
          {
            hydroVars.index = index;

            D3Q15::CalculateDensityVelocityFEq(hydroVars.f,
                                               hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z,
                                               hydroVars.f_eq);

            /*
             *  TODO optimise, at this point hydroVars.f_neq has not been computed yet, so we
             *  need to do it here for the shear-rate calculator. However, the streamer will do
             *  it again before DoCollide is called. Rearrange the code perhaps.
             */
            distribn_t f_neq[D3Q15::NUMVECTORS];
            for (unsigned f_index = 0; f_index < D3Q15::NUMVECTORS; f_index++)
            {
              f_neq[f_index] = hydroVars.f[f_index] - hydroVars.f_eq[f_index];
            }

            assert(hydroVars.index < (size_t) mTau.size());
            double old_tau_value = mTau[hydroVars.index];

            double shear_rate = D3Q15::CalculateShearRate(old_tau_value,
                                                          f_neq,
                                                          mTimeStep,
                                                          mSpaceStep,
                                                          hydroVars.density);

            mTau[hydroVars.index] = tNonNewtonianModel::CalculateTauForShearRate(shear_rate,
                                                                                 hydroVars.density,
                                                                                 mSpaceStep,
                                                                                 mTimeStep);
            // In some rheology models viscosity tends to infinity as shear rate goes to zero.
            assert( !std::isinf(mTau[hydroVars.index]) );
          }

          void DoCalculateFeq(HydroVars<LBGKNN>& hydroVars, site_t index)
          {
            D3Q15::CalculateFeq(hydroVars.density,
                                hydroVars.v_x,
                                hydroVars.v_y,
                                hydroVars.v_z,
                                hydroVars.f_eq);
          }

          distribn_t DoCollide(const LbmParameters* const lbmParams,
                               HydroVars<LBGKNN>& hydroVars,
                               unsigned int direction)
          {
            assert(hydroVars.index < (size_t) mTau.size());
            double omega = -1.0/mTau[hydroVars.index];
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
          /* Vector containing the current relaxation time for each site in the domain */
          std::vector<distribn_t> mTau;

          /* Current time step */
          distribn_t mTimeStep;

          /* Current space step */
          distribn_t mSpaceStep;
      };
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGKNN_H */
